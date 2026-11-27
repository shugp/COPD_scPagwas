"""Write a native editable Word letter, a plain coloc worksheet, and RIS references."""
import argparse, json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from xml.sax.saxutils import escape
import pandas as pd

REPO=Path(__file__).resolve().parents[1]

def reference_text(references):
    result=[]
    for i,r in enumerate(references,1):
        authors=[]
        for author in r["authors"][:3]:
            surname,given=author.split(", ",1)
            authors.append(surname+" "+"".join(w[0] for w in given.replace(".","").split()))
        names=", ".join(authors)+(", et al." if len(r["authors"])>3 else ".")
        journal=r.get("journal_abbreviation",r["journal"])
        result.append(f'[{i}] {names} {r["title"]}. {journal}. {r["year"]};{r["volume"]}({r["issue"]}):{r["pages"]}. https://doi.org/{r["doi"]}')
    return result

def write_ris(references,out):
    records=[]
    for i,r in enumerate(references,1):
        start,*end=r["pages"].split("–")
        fields=[("TY","JOUR"),("ID",f"MyeloidLetterRef{i:02}"),*[("AU",a) for a in r["authors"]],
                ("TI",r["title"]),("JO",r["journal"]),("PY",r["year"]),("VL",r["volume"]),
                ("IS",r["issue"]),("SP",start)]
        if end:fields.append(("EP",end[0]))
        fields.extend([("DO",r["doi"]),("UR","https://doi.org/"+r["doi"])])
        if r.get("pmid"):fields.append(("AN",r["pmid"]))
        records.append("\n".join(f"{key}  - {value}" for key,value in fields)+"\nER  -\n")
    (out/"Myeloid_COPD_References.ris").write_text("\n".join(records),encoding="utf-8")

def write_table(data_root,out):
    raw = pd.read_csv(data_root/"coloc.CSV")
    raw = raw[raw['PP.H4.abf'] > .6].sort_values(['Gene','cell_type','population','PP.H4.abf'],
                                             ascending=[True,True,True,False], kind='stable')
    headers = ['original_row_id'] + list(raw.columns[1:])
    def colname(i):
        s=''
        while i:
            i,r=divmod(i-1,26);s=chr(65+r)+s
        return s
    def cell_xml(v, address):
        if pd.isna(v): return f'<c r="{address}"/>'
        if isinstance(v, (int,float)):
            return f'<c r="{address}"><v>{v!r}</v></c>'
        return f'<c r="{address}" t="inlineStr"><is><t xml:space="preserve">{escape(str(v))}</t></is></c>'
    rows=[headers]+raw.values.tolist()
    sheet_rows=[]
    for i,row in enumerate(rows,1):
        sheet_rows.append(f'<row r="{i}">'+''.join(cell_xml(v,f'{colname(j)}{i}') for j,v in enumerate(row,1))+'</row>')
    ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    widths=[18,10,18,18,18,18,18,26,19,43,16,14]
    sheet=f'<worksheet xmlns="{ns}"><cols>'+''.join(f'<col min="{i}" max="{i}" width="{w}" customWidth="1"/>' for i,w in enumerate(widths,1))+'</cols><sheetData>'+''.join(sheet_rows)+'</sheetData></worksheet>'
    with ZipFile(out/'Myeloid_Coloc_Candidate_Pathways.xlsx','w',ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml','<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>')
        z.writestr('_rels/.rels','<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        z.writestr('xl/workbook.xml',f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="coloc_PP.H4_gt_0.6" sheetId="1" r:id="rId1"/></sheets></workbook>')
        z.writestr('xl/_rels/workbook.xml.rels','<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr('xl/worksheets/sheet1.xml',sheet)


def write_letter(manuscript,references,out,repository_url=None):
    title=manuscript["title"]
    body=manuscript["body"]
    legend=manuscript["legend"]
    refs=reference_text(references)
    availability=manuscript["availability"]
    code_link=repository_url or "［待作者上传GitHub后填写公开仓库URL；当前尚未公开］"
    availability=availability.replace("{repository_url}",code_link)
    W='http://schemas.openxmlformats.org/wordprocessingml/2006/main'
    def paragraph(text,style='Normal',bold=False):
     ppr=f'<w:pPr><w:pStyle w:val="{style}"/></w:pPr>'
     rpr='<w:rPr>'+('<w:b/>' if bold else '')+'</w:rPr>'
     return f'<w:p>{ppr}<w:r>{rpr}<w:t xml:space="preserve">{escape(text)}</w:t></w:r></w:p>'
    sections=[paragraph(title,'Title')]+[paragraph(t) for t in body]
    availability_paragraph=paragraph(availability,'Caption').replace('</w:pPr>','<w:keepLines/></w:pPr>')
    sections += [paragraph('数据与代码可获得性','Heading1'),availability_paragraph]
    sections += [paragraph('图注','Heading1')]+[paragraph(t,'Caption') for t in legend]
    sections += [paragraph('参考文献','Heading1')]+[paragraph(t,'Reference') for t in refs]
    section='<w:sectPr><w:footerReference w:type="default" r:id="rId2"/><w:pgSz w:w="11906" w:h="16838"/><w:pgMar w:top="1134" w:right="1134" w:bottom="1134" w:left="1134" w:header="500" w:footer="500"/></w:sectPr>'
    document=f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document xmlns:w="{W}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><w:body>'+''.join(sections)+section+'</w:body></w:document>'
    styles=f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
    <w:styles xmlns:w="{W}"><w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="Noto Serif CJK SC"/><w:sz w:val="22"/><w:lang w:val="en-US" w:eastAsia="zh-CN"/></w:rPr></w:rPrDefault><w:pPrDefault><w:pPr><w:spacing w:after="60" w:line="270" w:lineRule="auto"/><w:widowControl/></w:pPr></w:pPrDefault></w:docDefaults>
    <w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:pPr><w:jc w:val="both"/></w:pPr></w:style>
    <w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:basedOn w:val="Normal"/><w:pPr><w:keepNext/><w:jc w:val="center"/><w:spacing w:after="240"/></w:pPr><w:rPr><w:b/><w:rFonts w:eastAsia="Noto Sans CJK SC"/><w:sz w:val="32"/></w:rPr></w:style>
    <w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="Heading 1"/><w:basedOn w:val="Normal"/><w:pPr><w:keepNext/><w:spacing w:before="160" w:after="80"/></w:pPr><w:rPr><w:b/><w:sz w:val="24"/></w:rPr></w:style>
    <w:style w:type="paragraph" w:styleId="Caption"><w:name w:val="Caption"/><w:basedOn w:val="Normal"/><w:pPr><w:spacing w:line="264" w:after="80"/></w:pPr><w:rPr><w:sz w:val="20"/></w:rPr></w:style>
    <w:style w:type="paragraph" w:styleId="Reference"><w:name w:val="Reference"/><w:basedOn w:val="Normal"/><w:pPr><w:spacing w:line="220" w:after="40"/><w:jc w:val="left"/></w:pPr><w:rPr><w:sz w:val="19"/></w:rPr></w:style></w:styles>'''
    content_types='''<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/><Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/><Override PartName="/word/footer1.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/><Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/></Types>'''
    rels='''<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/></Relationships>'''
    docrels='''<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/><Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" Target="footer1.xml"/></Relationships>'''
    footer=f'''<?xml version="1.0" encoding="UTF-8"?><w:ftr xmlns:w="{W}"><w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:rPr><w:sz w:val="18"/></w:rPr><w:fldChar w:fldCharType="begin"/></w:r><w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r><w:r><w:fldChar w:fldCharType="end"/></w:r></w:p></w:ftr>'''
    core=f'''<?xml version="1.0" encoding="UTF-8"?><cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/"><dc:title>{escape(title)}</dc:title><dc:language>zh-CN</dc:language><dc:description>中文研究短函，包含正文、综合图图注及参考文献。</dc:description></cp:coreProperties>'''
    with ZipFile(out/'Myeloid_COPD_Letter_CN.docx','w',ZIP_DEFLATED) as z:
     for path,value in {'[Content_Types].xml':content_types,'_rels/.rels':rels,'word/document.xml':document,
                        'word/_rels/document.xml.rels':docrels,'word/styles.xml':styles,'word/footer1.xml':footer,
                        'docProps/core.xml':core}.items():z.writestr(path,value)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-root",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--repository-url")
    a=p.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    references=json.loads((REPO/"resources/references.json").read_text())
    manuscript=json.loads((REPO/"resources/manuscript_cn.json").read_text())
    write_table(a.data_root.resolve(),out)
    write_ris(references,out)
    write_letter(manuscript,references,out,a.repository_url)

if __name__=="__main__":main()
