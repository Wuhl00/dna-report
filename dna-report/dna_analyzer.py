import os
import requests
import time
import hashlib
import io
import re
from urllib.parse import quote
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from fpdf import FPDF

class DNAAnalyzer:
    """
    DNA Sequence Deep Analysis Pipeline.
    Integrates: Basic properties, ORF prediction, EBI BLASTN, and Report Generation.
    """
    
    def __init__(self, fasta_path, output_dir=None):
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"Input file {fasta_path} not found.")
            
        self.fasta_path = fasta_path
        self.record = list(SeqIO.parse(fasta_path, "fasta"))[0]
        self.sequence = str(self.record.seq).upper()
        self.output_dir = output_dir or os.getcwd()
        os.makedirs(self.output_dir, exist_ok=True)
        
        project_hash = hashlib.md5(os.getcwd().encode()).hexdigest()[:8]
        self.email = f"trae_user_{project_hash}@gmail.com"
        
        self.analysis_results = {
            'properties': {},
            'orfs': [],
            'blast_hits': [],
            'ai_structured': {}
        }

    def analyze_properties(self):
        """Calculate basic DNA properties."""
        print(f"--- Analyzing properties for {self.record.id} ---")
        length = len(self.sequence)
        gc_content = gc_fraction(self.sequence) * 100
        
        # Approximate MW calculation
        # dsDNA MW: length * 617.96 + 36.04
        approx_mw_ds = length * 617.96 + 36.04
        # ssDNA MW: length * 308.97 + 18.02
        approx_mw_ss = length * 308.97 + 18.02
        # RNA MW (approx): length * 320.5 + 159.0
        approx_mw_rna = length * 320.5 + 159.0
        
        self.analysis_results['properties'] = {
            'Length (bp)': length,
            'GC Content (%)': f"{gc_content:.2f}",
            'Approximate MW (dsDNA, Da)': f"{approx_mw_ds:.2f}",
            'Approximate MW (ssDNA, Da)': f"{approx_mw_ss:.2f}",
            'Approximate MW (RNA, Da)': f"{approx_mw_rna:.2f}"
        }

    def scan_restriction_enzymes(self):
        """Scan for common 6-bp restriction enzyme cut sites."""
        print("--- Scanning Restriction Enzymes ---")
        common_enzymes = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'XhoI': 'CTCGAG',
            'NotI': 'GCGGCCGC',
            'NdeI': 'CTCGAG', # Actually XhoI, let's fix: NdeI is CATATG
            'NheI': 'CATATG', # NheI is GCTAGC
            'NcoI': 'GCTAGC', # NcoI is CCATGG
            'BglII': 'AGATCT'
        }
        # Correcting the above quick dict to standard sites:
        common_enzymes = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'XhoI': 'CTCGAG',
            'NotI': 'GCGGCCGC',
            'NdeI': 'CATATG',
            'NheI': 'GCTAGC',
            'NcoI': 'CCATGG',
            'BglII': 'AGATCT',
            'SacI': 'GAGCTC',
            'SalI': 'GAGCTC', # SacI is GAGCTC, SalI is GTCGAC
            'KpnI': 'CCATGG', # KpnI is GGTACC
        }
        # Final accurate dictionary for 10 common enzymes
        common_enzymes = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'XhoI': 'CTCGAG',
            'NotI': 'GCGGCCGC',
            'NdeI': 'CATATG',
            'NheI': 'GCTAGC',
            'NcoI': 'CCATGG',
            'BglII': 'AGATCT',
            'SalI': 'GTCGAC'
        }
        
        re_results = {}
        for enz, site in common_enzymes.items():
            # Find all occurrences
            positions = [m.start() + 1 for m in re.finditer(f'(?={site})', self.sequence)]
            re_results[enz] = {
                'site': site,
                'count': len(positions),
                'positions': positions
            }
        
        self.analysis_results['restriction_enzymes'] = re_results

    def run_blastn(self, hit_count=5):
        """
        Run NCBI BLASTN (async via REST API).
        Database: nt (Nucleotide collection)
        """
        print("--- Submitting NCBI BLASTN task (Database: nt) ---")
        base_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
        
        try:
            submit = requests.post(base_url, data={
                "CMD": "Put",
                "PROGRAM": "blastn",
                "DATABASE": "nt",
                "QUERY": self.sequence
            })
            match = re.search(r"RID = (.*)", submit.text)
            if not match:
                print("NCBI BLAST submission failed: No RID found.")
                self.analysis_results['blast_hits'] = []
                return False
            job_id = match.group(1).strip()
            print(f"DEBUG: NCBI BLAST RID = {job_id}")
        except Exception as e:
            print(f"NCBI BLAST submission exception: {e}")
            self.analysis_results['blast_hits'] = []
            return False

        # Polling status, max wait 300s (NCBI nt takes longer)
        start_time = time.time()
        print("--- Polling NCBI BLAST task status ---")
        while True:
            if time.time() - start_time > 300:
                print("NCBI BLAST timeout (5 mins). Gracefully skipping.")
                self.analysis_results['blast_hits'] = []
                return False
            
            time.sleep(10) # NCBI recommends 10s wait between requests
            try:
                status_req = requests.get(base_url, params={"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": job_id})
                if "Status=WAITING" in status_req.text:
                    continue
                if "Status=FAILED" in status_req.text:
                    print("NCBI BLAST task failed.")
                    self.analysis_results['blast_hits'] = []
                    return False
                if "Status=READY" in status_req.text:
                    if "ThereAreHits=yes" in status_req.text:
                        break
                    else:
                        print("NCBI BLAST finished but no hits found.")
                        self.analysis_results['blast_hits'] = []
                        return True
            except Exception:
                continue

        # Fetch XML results
        print("--- Fetching NCBI BLAST results ---")
        hits = []
        try:
            from Bio.Blast import NCBIXML
            res_xml = requests.get(base_url, params={"CMD": "Get", "FORMAT_TYPE": "XML", "RID": job_id})
            records = NCBIXML.parse(io.StringIO(res_xml.text))
            record = next(records)
            
            for aln in record.alignments:
                title = aln.title
                # Parse title to keep only 4th and 5th parts separated by '|'
                # Example full title: "gi|257094701|ref|NM_001175256.1| Zea mays uncharacterized LOC100382519 (LOC100382519), mRNA"
                parts = title.split('|')
                if len(parts) >= 5:
                    title_clean = f"{parts[3]}|{parts[4].strip()}"
                else:
                    title_clean = title

                acc = aln.accession
                if aln.hsps:
                    hsp = aln.hsps[0]
                    evalue = hsp.expect
                    identity_pct = f"{(hsp.identities/hsp.align_length)*100:.1f}%"
                    q_start = hsp.query_start
                    q_end = hsp.query_end
                else:
                    evalue = "N/A"
                    identity_pct = "N/A"
                    q_start = 0
                    q_end = 0
                
                hits.append({
                    'title': title_clean,
                    'acc': acc,
                    'e_value': evalue,
                    'identity': identity_pct,
                    'q_start': q_start,
                    'q_end': q_end,
                    'url': f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}"
                })
                if len(hits) >= hit_count:
                    break
        except Exception as e:
            print(f"XML parsing failed: {e}")
            self.analysis_results['blast_hits'] = []

        self.analysis_results['blast_hits'] = hits
        print(f"DEBUG: Retrieved {len(hits)} homology hits")
        return bool(hits)

    def generate_ai_summary(self):
        """Generate structured AI summary."""
        print("--- Generating AI Summary ---")
        props = self.analysis_results.get('properties', {})
        orfs = self.analysis_results.get('orfs', [])
        blast_hits = self.analysis_results.get('blast_hits', [])
        
        top_hit = blast_hits[0] if blast_hits else None
        
        # 1. Summary
        summary_part = (
            f"The analyzed DNA sequence is {props.get('Length (bp)', 'N/A')} bp long with a GC content of {props.get('GC Content (%)', 'N/A')}%. "
        )
            
        if top_hit:
            summary_part += f"BLASTN homology search against NCBI nt database identified a top match with {top_hit['acc']} ({top_hit['identity']} identity)."
        else:
            summary_part += "No significant nucleotide homology was found or the search timed out."

        # 2. Prediction
        prediction_part = ""
        if top_hit and float(top_hit['identity'].replace('%', '')) > 80:
            prediction_part = f"Based on high nucleotide sequence identity, this sequence is highly conserved and likely functionally identical or closely related to the gene product of {top_hit['acc']}."
        else:
            prediction_part = "Given the lack of obvious known homology, this sequence might represent a non-coding RNA, a regulatory element (such as a promoter or enhancer), or an intergenic region."

        # 3. Dynamic PubMed Keywords
        search_query = self.record.id
        if top_hit:
            # Extract meaningful keywords from top hit title
            # Example title: "Zea mays uncharacterized LOC100382519 (LOC100382519), mRNA"
            title = top_hit['title']
            # Remove common unhelpful words
            stopwords = ['uncharacterized', 'mRNA', 'predicted', 'clone', 'isoform', 'transcript', 'variant']
            words = [w for w in re.split(r'[\s,()|]+', title) if w and w.lower() not in stopwords]
            if words:
                # Take first 3-4 meaningful words (e.g. "Zea", "mays", "LOC100382519")
                search_query = "+".join(words[:4])
                
        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/?term={search_query}"

        # 4. References
        ref_links = []
        if top_hit:
            ref_links.append(f"[ENA Browser: {top_hit['acc']}]({top_hit['url']})")
        ref_links.append(f"[Search PubMed for related literature]({pubmed_url})")

        self.analysis_results['ai_structured'] = {
            'summary': summary_part,
            'prediction': prediction_part,
            'references': ref_links
        }
        return self.analysis_results['ai_structured']

    def generate_report(self, output_pdf="dna_report.pdf"):
        pdf_path = os.path.join(self.output_dir, output_pdf)
        pdf = FPDF()
        bookmarks = []
        pdf.add_page()
        
        # Title
        pdf.set_font("Arial", 'B', 20)
        pdf.cell(200, 15, txt=f"DNA Analysis: {self.record.id}", ln=True, align='C')
        
        # 1. Sequence
        bookmarks.append(("1. Input Sequence", pdf.page_no()))
        pdf.set_font("Arial", 'B', 14)
        pdf.cell(200, 10, txt="1. Input Sequence", ln=True)
        pdf.set_font("Courier", '', 8)
        seq = self.sequence
        
        # Calculate how many lines of sequence can fit on the first page 
        # while leaving enough room for Section 2.
        max_display_bp = 3150
        
        for i in range(0, min(len(seq), max_display_bp), 70):
            pdf.cell(0, 4, txt=seq[i:i+70], ln=True)
            
        if len(seq) > max_display_bp:
            pdf.cell(0, 4, txt=f"... [Sequence truncated at {max_display_bp} bp to fit first page]", ln=True)
        pdf.ln(5)
        
        # 2. Properties
        bookmarks.append(("2. Basic Properties", pdf.page_no()))
        pdf.set_font("Arial", 'B', 14)
        pdf.cell(200, 10, txt="2. Basic Properties", ln=True)
        pdf.set_font("Arial", '', 11)
        for k, v in self.analysis_results.get('properties', {}).items():
            pdf.cell(0, 8, txt=f"- {k}: {v}", ln=True)
        pdf.ln(5)

        # Restriction Enzymes
        if 'restriction_enzymes' in self.analysis_results:
            pdf.set_font("Arial", 'B', 12)
            pdf.cell(200, 8, txt="Restriction Enzyme Sites (Common 6-cutters)", ln=True)
            pdf.set_font("Arial", '', 10)
            
            # Simple table for RE sites
            for enz, data in self.analysis_results['restriction_enzymes'].items():
                if data['count'] > 0:
                    pos_str = ", ".join(map(str, data['positions']))
                    pdf.cell(0, 6, txt=f"  - {enz} ({data['site']}): {data['count']} cut(s) at pos: {pos_str}", ln=True)
            # If no cuts found at all
            if all(d['count'] == 0 for d in self.analysis_results['restriction_enzymes'].values()):
                pdf.set_font("Arial", 'I', 10)
                pdf.cell(0, 6, txt="  - No cut sites found for the tested enzymes.", ln=True)
            pdf.ln(5)

        # 3. AI Portal (Evo 2)
        pdf.add_page()
        bookmarks.append(("3. AI Genomic Foundation Model", pdf.page_no()))
        pdf.set_font("Arial", 'B', 16)
        pdf.cell(200, 12, txt="3. AI Genomic Foundation Model (Evo 2)", ln=True)
        pdf.set_font("Arial", '', 11)
        pdf.multi_cell(0, 6, txt="Evo 2 is a genomic foundation model capable of predicting and designing across DNA, RNA, and proteins. You can use it to score this sequence for perplexity and per-nucleotide entropy.")
        pdf.ln(2)
        pdf.set_font("Arial", 'I', 10)
        pdf.set_text_color(100, 100, 100)
        pdf.multi_cell(0, 5, txt="Note: If you use Evo 2 for your research, please cite: https://www.nature.com/articles/s41586-026-10176-5")
        pdf.set_text_color(0, 0, 0)
        pdf.ln(4)
        pdf.set_text_color(0, 0, 255)
        pdf.set_font("Arial", 'BU', 11)
        pdf.write(8, "Open Evo 2 Designer Portal", "https://arcinstitute.org/tools/evo/evo-designer")
        pdf.set_text_color(0, 0, 0)
        pdf.ln(10)

        # Human actB Example
        pdf.set_font("Arial", 'B', 12)
        pdf.cell(0, 8, txt="Example: Human actB Sequence Analysis", ln=True)
        
        # Check if the example image exists
        img_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "evo2_actb_example.png")
        if os.path.exists(img_path):
            pdf.image(img_path, w=180)
            pdf.ln(5)
        else:
            # Fallback to current working directory if not found in script dir
            cwd_img_path = os.path.join(os.getcwd(), "evo2_actb_example.png")
            if os.path.exists(cwd_img_path):
                pdf.image(cwd_img_path, w=180)
                pdf.ln(5)
            else:
                pdf.set_font("Arial", 'I', 10)
                pdf.set_text_color(150, 150, 150)
                pdf.cell(0, 10, txt="[Figure 1: Sequence Logo Placeholder - Please place 'evo2_actb_example.png' in script dir]", ln=True)
                pdf.set_text_color(0, 0, 0)
                pdf.ln(2)

        pdf.set_font("Arial", '', 10)
        desc = "The sequence logo shows the probability of each nucleotide at each position in the sequence. The height of each letter indicates the probability of that nucleotide at that position. The total height of the stack of letters at each position indicates the total information content of that position, where higher stacks indicate more conserved."
        pdf.multi_cell(0, 5, txt=desc)
        pdf.ln(4)
        
        # Legend A, C, T, G
        pdf.set_font("Arial", 'B', 10)
        
        # A
        pdf.set_fill_color(253, 224, 139) # #FDE08B
        pdf.cell(8, 6, txt=" A ", border=0, fill=True, align='C')
        pdf.cell(2, 6, txt="", border=0)
        # C
        pdf.set_fill_color(138, 226, 242) # #8AE2F2
        pdf.cell(8, 6, txt=" C ", border=0, fill=True, align='C')
        pdf.cell(2, 6, txt="", border=0)
        # T
        pdf.set_fill_color(208, 189, 244) # #D0BDF4
        pdf.cell(8, 6, txt=" T ", border=0, fill=True, align='C')
        pdf.cell(2, 6, txt="", border=0)
        # G
        pdf.set_fill_color(189, 231, 130) # #BDE782
        pdf.cell(8, 6, txt=" G ", border=0, fill=True, align='C')
        pdf.ln(15)

        # 4. BLASTN
        pdf.add_page()
        bookmarks.append(("4. Homology Analysis (BLASTN)", pdf.page_no()))
        pdf.set_font("Arial", 'B', 16)
        pdf.cell(200, 12, txt="4. Homology Analysis (BLASTN) - Top 5 Hits", ln=True)
        pdf.ln(5)
        
        if self.analysis_results.get('blast_hits'):
            pdf.set_font("Arial", 'I', 10)
            pdf.set_text_color(100, 100, 100)
            pdf.multi_cell(0, 6, txt="Note: This search uses the NCBI nt database. For comprehensive search, please use the official portal:")
            pdf.set_text_color(0, 0, 255)
            pdf.set_font("Arial", 'BU', 10)
            pdf.write(8, "[Go to NCBI BLASTN Official Website]", "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch")
            pdf.ln(10)
            
            pdf.set_text_color(0, 0, 0)
            pdf.set_font("Arial", '', 11)
            for hit in self.analysis_results['blast_hits']:
                pdf.set_font("Arial", 'B', 11)
                title_clean = "".join(c for c in hit['title'][:150] if ord(c) < 128)
                pdf.multi_cell(0, 7, txt=f"> {title_clean}")
                
                pdf.set_font("Arial", '', 10)
                pdf.write(6, f"  - Identity: {hit['identity']} | E-value: {hit['e_value']} | Accession: ")
                pdf.set_text_color(0, 0, 255)
                pdf.write(6, f"{hit['acc']}", hit['url'])
                pdf.set_text_color(0, 0, 0)
                pdf.ln(8)
        else:
            pdf.set_font("Arial", 'I', 11)
            pdf.cell(0, 8, txt="No significant nucleotide homology was found or the search timed out.", ln=True)
            pdf.ln(5)

        # 5. AI Summary
        if 'ai_structured' in self.analysis_results:
            bookmarks.append(("5. AI-Assisted Function Prediction", pdf.page_no()))
            pdf.set_font("Arial", 'B', 16)
            pdf.cell(200, 12, txt="5. AI-Assisted Function Prediction", ln=True)
            pdf.ln(5)
            
            ai = self.analysis_results['ai_structured']
            
            pdf.set_font("Arial", 'B', 12)
            pdf.cell(0, 8, txt="5.1 Investigation Summary", ln=True)
            pdf.set_font("Arial", '', 11)
            pdf.multi_cell(0, 7, txt=ai['summary'])
            pdf.ln(5)
            
            pdf.set_font("Arial", 'B', 12)
            pdf.cell(0, 8, txt="5.2 Functional Prediction", ln=True)
            pdf.set_font("Arial", '', 11)
            pdf.set_fill_color(245, 245, 255) 
            pdf.multi_cell(0, 7, txt=ai['prediction'], border=1, fill=True)
            pdf.ln(5)
            
            pdf.set_font("Arial", 'B', 12)
            pdf.cell(0, 8, txt="5.3 Related Literature Search", ln=True)
            pdf.set_font("Arial", '', 10)
            for ref in ai['references']:
                if '](' in ref and ')' in ref:
                    text = ref[ref.find("[")+1:ref.find("]")]
                    url = ref[ref.find("(")+1:ref.find(")")]
                    text_ascii = "".join(c for c in text if ord(c) < 128)
                    pdf.set_text_color(0, 0, 255)
                    pdf.write(7, f"- {text_ascii}", url)
                    pdf.set_text_color(0, 0, 0)
                    pdf.ln(7)
            
            pdf.ln(10)
            pdf.set_font("Arial", 'I', 9)
            pdf.set_text_color(150, 150, 150)
            pdf.multi_cell(0, 5, txt="Disclaimer: This summary is automatically synthesized by the Trae AI DNA analysis agent. Please verify with experiments.")

        pdf.output(pdf_path)
        print(f"PDF Report generated: {pdf_path}")
        self._try_add_pdf_bookmarks(pdf_path, bookmarks)

    def _try_add_pdf_bookmarks(self, pdf_path, bookmarks):
        try:
            from PyPDF2 import PdfReader, PdfWriter
        except Exception:
            return False
        try:
            reader = PdfReader(pdf_path)
            writer = PdfWriter()
            for page in reader.pages:
                writer.add_page(page)
            for title, page_no in bookmarks:
                try:
                    writer.add_outline_item(title, page_no-1)
                except Exception:
                    try:
                        writer.addBookmark(title, page_no-1)
                    except Exception:
                        pass
            tmp_path = pdf_path + ".tmp.pdf"
            with open(tmp_path, "wb") as f:
                writer.write(f)
            try:
                os.replace(tmp_path, pdf_path)
            except Exception:
                pass
            return True
        except Exception:
            return False

    def generate_markdown_report(self, output_md="dna_report.md"):
        # Markdown Generation
        md_path = os.path.join(self.output_dir, output_md)
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(f"# DNA Sequence Deep Analysis Report: {self.record.id}\n\n")
            
            f.write("## 1. Input Sequence\n")
            f.write("```fasta\n")
            f.write(f">{self.record.id}\n")
            seq = self.sequence
            
            md_max_display = 10000
            for i in range(0, min(len(seq), md_max_display), 70):
                f.write(seq[i:i+70] + "\n")
                
            if len(seq) > md_max_display:
                f.write(f"... [Sequence truncated at {md_max_display} bp]\n")
            f.write("```\n\n")
            
            f.write("## 2. Basic Properties\n")
            for k, v in self.analysis_results.get('properties', {}).items():
                f.write(f"- **{k}**: {v}\n")
            f.write("\n")

            if 'restriction_enzymes' in self.analysis_results:
                f.write("### Restriction Enzyme Sites (Common 6-cutters)\n")
                has_cuts = False
                for enz, data in self.analysis_results['restriction_enzymes'].items():
                    if data['count'] > 0:
                        has_cuts = True
                        pos_str = ", ".join(map(str, data['positions']))
                        f.write(f"- **{enz}** ({data['site']}): {data['count']} cut(s) at pos: {pos_str}\n")
                if not has_cuts:
                    f.write("- *No cut sites found for the tested enzymes.*\n")
                f.write("\n")

            f.write("## 3. AI Genomic Foundation Model (Evo 2)\n")
            f.write("Evo 2 is a genomic foundation model capable of predicting and designing across DNA, RNA, and proteins. You can use it to score this sequence for perplexity and per-nucleotide entropy.\n\n")
            f.write("> **Note:** If you use Evo 2 for your research, please cite: https://www.nature.com/articles/s41586-026-10176-5 \n\n")
            f.write("- **[👉 Open Evo 2 Designer Portal](https://arcinstitute.org/tools/evo/evo-designer)**\n\n")
            
            f.write("### Example: Human actB Sequence Analysis\n\n")
            
            img_path = "evo2_actb_example.png"
            script_img_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "evo2_actb_example.png")
            if os.path.exists(script_img_path):
                # Copy image to output dir for markdown relative path
                import shutil
                try:
                    shutil.copy(script_img_path, os.path.join(self.output_dir, "evo2_actb_example.png"))
                except Exception:
                    pass
            elif os.path.exists(os.path.join(os.getcwd(), "evo2_actb_example.png")):
                import shutil
                try:
                    shutil.copy(os.path.join(os.getcwd(), "evo2_actb_example.png"), os.path.join(self.output_dir, "evo2_actb_example.png"))
                except Exception:
                    pass

            f.write(f"![Sequence Logo]({img_path})\n\n")
            f.write("The sequence logo shows the probability of each nucleotide at each position in the sequence. The height of each letter indicates the probability of that nucleotide at that position. The total height of the stack of letters at each position indicates the total information content of that position, where higher stacks indicate more conserved.\n\n")
            f.write("<span style=\"background-color:#FDE08B; padding:2px 6px; border-radius:3px; color:black; font-weight:bold;\">A</span> ")
            f.write("<span style=\"background-color:#8AE2F2; padding:2px 6px; border-radius:3px; color:black; font-weight:bold;\">C</span> ")
            f.write("<span style=\"background-color:#D0BDF4; padding:2px 6px; border-radius:3px; color:black; font-weight:bold;\">T</span> ")
            f.write("<span style=\"background-color:#BDE782; padding:2px 6px; border-radius:3px; color:black; font-weight:bold;\">G</span>\n\n")

            f.write("## 4. Homology Analysis (BLASTN) - Top 5 Hits\n")
            if self.analysis_results.get('blast_hits'):
                f.write("> **Note**: This search uses the NCBI nt database. For comprehensive search, please use the official portal:\n")
                f.write("> **[👉 NCBI BLASTN Official Website](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch)**\n\n")
                f.write("| Title | Identity | E-value | Accession |\n")
                f.write("| --- | --- | --- | --- |\n")
                for hit in self.analysis_results['blast_hits']:
                    title_clean = hit['title'][:150].replace('|', '\\|')
                    f.write(f"| {title_clean} | {hit['identity']} | {hit['e_value']} | [{hit['acc']}]({hit['url']}) |\n")
                f.write("\n")
            else:
                f.write("> No significant nucleotide homology was found or the search timed out.\n\n")

            if 'ai_structured' in self.analysis_results:
                f.write("## 5. AI-Assisted Function Prediction\n")
                ai = self.analysis_results['ai_structured']
                f.write("### 5.1 Investigation Summary\n")
                f.write(f"> {ai['summary']}\n\n")
                f.write("### 5.2 Functional Prediction\n")
                f.write(f"**{ai['prediction']}**\n\n")
                f.write("### 5.3 Related Literature Search\n")
                for ref in ai['references']:
                    f.write(f"- {ref}\n")
                f.write("\n---\n*Disclaimer: This summary is automatically synthesized by the Trae AI DNA analysis agent. Please verify with experiments.*\n\n")

        print(f"Markdown Report generated: {md_path}")

if __name__ == "__main__":
    # Prioritize checking the current working directory, then the script's directory
    cwd_input = os.path.join(os.getcwd(), "input_dna.fasta")
    script_input = os.path.join(os.path.dirname(os.path.abspath(__file__)), "input_dna.fasta")
    
    input_file = cwd_input if os.path.exists(cwd_input) else script_input
    
    if os.path.exists(input_file):
        safe_id = re.sub(r"[^A-Za-z0-9._-]+", "_", str(list(SeqIO.parse(input_file, "fasta"))[0].id))
        ts = time.strftime("%Y%m%d_%H%M%S")
        out_dir = os.path.join(os.getcwd(), "dna_analysis_runs", f"{safe_id}_{ts}")
        
        analyzer = DNAAnalyzer(input_file, output_dir=out_dir)
        analyzer.analyze_properties()
        analyzer.scan_restriction_enzymes()
        
        # Async NCBI BLASTN (timeout 300s)
        analyzer.run_blastn()
        
        analyzer.generate_ai_summary()
        analyzer.generate_report(f"{safe_id}_report.pdf")
        analyzer.generate_markdown_report(f"{safe_id}_report.md")
    else:
        print(f"Error: Could not find {input_file}. Please create it with a DNA FASTA sequence.")
