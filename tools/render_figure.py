"""Regenerate manuscript Figure 1 from existing numerical exports, without fitting models."""
import argparse
from pathlib import Path
import os, tempfile, subprocess, copy
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "myeloid-mpl"))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import xml.etree.ElementTree as ET

def pathway_context(base):
    """Keep module enrichment and Reactome membership distinct; no coloc probability plot."""
    source=base/"disease_pathways"
    modules=pd.read_csv(source/"candidate_module_reactome_enrichment.csv")
    sets={}
    for line in (source/"Reactome_2022.gmt").read_text().splitlines():
        fields=line.split("\t");sets[fields[0]]=set(fields[2:])
    candidates=pd.read_csv(base/"results_cleaned/coloc_overall_gene_tests.csv")["Gene"].dropna().unique()
    rows=[]
    for gene in candidates:
        for pid,label in [("R-HSA-983170","MHC-I肽装载"),("R-HSA-2132295","MHC-II抗原呈递"),
                          ("R-HSA-1236975","交叉呈递"),("R-HSA-877300","IFN-γ信号")]:
            term=next(p for p in sets if pid in p)
            hit=modules[modules.candidate.eq(gene)&modules.pathway.eq(term)]
            symbols={gene,"LRMDA"} if gene=="C10orf11" else {gene}
            rows.append({"候选基因":gene,"通路":label,"共表达模块富集FDR":hit.FDR.iloc[0] if len(hit) else np.nan,
                         "Reactome2022成员注释":"已收录" if symbols&sets[term] else "未收录"})
    return pd.DataFrame(rows)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    root=args.data_root.resolve();base=root/"Myeloid";out=args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="myeloid-panel-a-") as directory:
        pdf_svg=Path(directory)/"original.svg"
        subprocess.run(["pdftocairo","-svg",str(root/"gPAS/scPagwas.inte_all.bar_positie_nagtive1.pdf"),str(pdf_svg)],check=True)
        plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.2,'axes.labelsize':7.2,'axes.titlesize':8,
         'xtick.labelsize':6.6,'ytick.labelsize':6.6,'axes.spines.top':False,'axes.spines.right':False,
         'svg.fonttype':'none','svg.hashsalt':'myeloid-letter-revision','lines.linewidth':.9})
        colors={'Healthy':'#4575A5','COPD':'#B65B4C'}
        fig=plt.figure(figsize=(14,15.6))
        gs=GridSpec(4,12,figure=fig,height_ratios=[3.15,2.9,3.25,2.7],left=.065,right=.965,bottom=.055,top=.97,hspace=.70,wspace=1.7)
        def panel(ax,letter,title):
         ax.set_title(title,loc='left',pad=6)
         ax.text(-.13,1.10,letter,transform=ax.transAxes,fontweight='bold',fontsize=11,va='top')
        # A reserves the exact source-PDF vector panel for direct SVG embedding.
        axA=fig.add_subplot(gs[0,:3]);axA.set_axis_off();panel(axA,'A','scPagwas positive / negative cells')
        # B. Original UMAP coordinates and five broad annotations.
        meta=pd.read_csv(base/'trajectory/cell_pseudotime_state.csv').set_index('cell')
        xy=pd.read_csv(base/'results_cleaned/umap_coordinates.csv',index_col=0).loc[meta.index].to_numpy()
        ax=fig.add_subplot(gs[0,3:6]);order=['Macrophage','Monocyte','DC','Mast/basophil','Neutrophil']
        for i,t in enumerate(order):
         mask=meta.myeloid_subtype.eq(t).to_numpy();ax.scatter(xy[mask,0],xy[mask,1],s=.18,color=plt.get_cmap('tab10')(i),linewidths=0,label=t)
        ax.set(xlabel='UMAP 1',ylabel='UMAP 2',aspect='equal');ax.legend(frameon=False,fontsize=5.8,markerscale=6,
         loc='upper left',bbox_to_anchor=(-.12,-.18),ncol=2,handlelength=.7,columnspacing=.6)
        panel(ax,'B','Myeloid subtypes')
        # C. Myeloid lineage markers.
        ax=fig.add_subplot(gs[0,6:]);dot=pd.read_csv(base/'results_cleaned/annotation_marker_expression.csv')
        genes=['CD68','C1QC','FABP4','FCN1','VCAN','FCGR3A','CD1C','CLEC4C','TCF4','TPSAB1','CPA3','HDC','FCGR3B','CSF3R','CXCR2']
        dot=dot[dot.gene.isin(genes)];rows=order[::-1]
        pts=ax.scatter(dot.gene.map({g:i for i,g in enumerate(genes)}),dot.subtype.map({g:i for i,g in enumerate(rows)}),
         s=dot.percent_expressing*.48,c=dot.mean_log1p_expression,cmap='Blues',edgecolors='none')
        ax.set_xticks(range(len(genes)),genes,rotation=60,ha='right');ax.set_yticks(range(5),rows);ax.set_xlim(-.6,len(genes)-.4)
        cb=fig.colorbar(pts,ax=ax,fraction=.025,pad=.03);cb.set_label('Mean log1p expression',fontsize=6.5)
        for pct in [25,50,100]:ax.scatter([],[],s=pct*.48,color='gray',label=str(pct)+'%')
        ax.legend(frameon=False,fontsize=6,loc='upper left',bbox_to_anchor=(0,-.32),ncol=3)
        panel(ax,'C','Representative annotation markers')
        # D. Matched disease-group trajectory views.
        a=meta[meta.trajectory_scope.eq('monocyte_macrophage')];u=xy[meta.trajectory_scope.eq('monocyte_macrophage')]
        sub=gs[1,:5].subgridspec(1,2,wspace=.2)
        for j,status in enumerate(['Healthy','COPD']):
         ax=fig.add_subplot(sub[0,j]);mask=a.status.eq(status).to_numpy()
         ax.scatter(u[mask,0],u[mask,1],c=a.pseudotime.to_numpy()[mask],cmap='viridis',vmin=0,vmax=1,s=.2,linewidths=0)
         ax.set(xlabel='UMAP 1',ylabel='UMAP 2' if j==0 else '',aspect='equal',xlim=(-9,13),ylim=(-10,14))
         ax.set_title(status,fontsize=7.5)
         if j==0:ax.text(-.16,1.12,'D',transform=ax.transAxes,fontweight='bold',fontsize=11)
         if j==1:
          cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap='viridis'),ax=ax,fraction=.045,pad=.035);cb.ax.set_title('DPT',fontsize=6.3,pad=4)
        # E. Previously generated GSEA results, shown as a disease-by-State enrichment matrix.
        ax=fig.add_subplot(gs[1,6:]);gsea=pd.read_csv(base/'disease_pathways/gsea_all_results.csv')
        defs={'MHC-I peptide loading':'Antigen Presentation: Folding, Assembly, Peptide Loading Of Class I MHC R-HSA-983170',
         'MHC-II presentation':'MHC Class II Antigen Presentation R-HSA-2132295',
         'Cross-presentation':'Antigen processing-Cross Presentation R-HSA-1236975',
         'IFN-gamma signaling':'Interferon Gamma Signaling R-HSA-877300','IFN-gamma response':'Interferon Gamma Response'}
        parts=[]
        for name,term in defs.items():
         d=gsea[gsea.Term.str.split(' | ',regex=False).str[-1].eq(term)&gsea.contrast.str.match(r'COPD_State[1-6]$')].copy()
         d['path']=name;d['State']=d.contrast.str.replace('COPD_','',regex=False);parts.append(d)
        g=pd.concat(parts);m=g.pivot(index='path',columns='State',values='NES').reindex(defs)
        q=g.pivot(index='path',columns='State',values='FDR q-val').reindex_like(m)
        im=ax.pcolormesh(np.arange(7)-.5,np.arange(6)-.5,m.to_numpy(),cmap='RdBu_r',vmin=-3,vmax=3,shading='flat');ax.invert_yaxis()
        for i in range(len(m)):
         for j in range(6):
          if q.iloc[i,j]<.05:ax.text(j,i,'*',ha='center',va='center',fontsize=9)
        ax.set_xticks(range(6),m.columns);ax.set_yticks(range(len(m)),m.index)
        cb=fig.colorbar(im,ax=ax,fraction=.025,pad=.03);cb.set_label('GSEA NES',fontsize=6.5)
        panel(ax,'E','COPD vs Healthy GSEA (* FDR < 0.05)')
        # F. Previous coexpression panel retained.
        ax=fig.add_subplot(gs[2,:5]);cor=pd.read_csv(base/'disease_pathways/antigen_residual_coexpression.csv')
        cands=['HLA-DQA2','HLA-DQB1','HLA-DQB2','HLA-DRB1','B2M','ERAP2']
        partners=['CD74','HLA-DMA','HLA-DMB','CTSS','TAP1','TAPBP','HLA-A','STAT1']
        m=cor.pivot(index='candidate',columns='partner',values='partial_correlation').loc[cands,partners]
        q=cor.pivot(index='candidate',columns='partner',values='FDR').loc[cands,partners]
        im=ax.pcolormesh(np.arange(9)-.5,np.arange(7)-.5,m.to_numpy(),cmap='RdBu_r',vmin=-1,vmax=1,shading='flat');ax.invert_yaxis()
        for i in range(6):
         for j in range(8):
          if q.iloc[i,j]<.05:ax.text(j,i,'*',ha='center',va='center',fontsize=8)
        ax.set_xticks(range(8),partners,rotation=50,ha='right');ax.set_yticks(range(6),cands)
        cb=fig.colorbar(im,ax=ax,fraction=.025,pad=.03);cb.set_label('Partial r',fontsize=6.7)
        panel(ax,'F','Coexpression adjusted for donor and State')
        # G. Candidate coexpression modules and antigen-presentation pathway context.
        aff=pathway_context(base)
        labels=['HLA-DQA2','HLA-DQB1','HLA-DQB2','HLA-DRB1','B2M','ERAP2','HSPA5','MAP1LC3B',
                'C10orf11','TUFM','LRPAP1','RP11-345J4.5','HEXA','SULT1A1','RBM6']
        paths=['MHC-I肽装载','MHC-II抗原呈递','交叉呈递','IFN-γ信号']
        ax=fig.add_subplot(gs[2,7:])
        mat=aff.pivot(index='候选基因',columns='通路',values='共表达模块富集FDR').loc[labels,paths]
        val=-np.log10(mat.clip(lower=1e-15));val=val.where(mat<.05,0).clip(upper=6)
        im=ax.pcolormesh(np.arange(5)-.5,np.arange(len(labels)+1)-.5,val.to_numpy(),cmap='Blues',vmin=0,vmax=6,shading='flat');ax.invert_yaxis()
        member=aff.pivot(index='候选基因',columns='通路',values='Reactome2022成员注释').loc[labels,paths]
        for i in range(len(labels)):
         for j in range(4):
          if member.iloc[i,j]=='已收录':ax.scatter(j,i,s=20,marker='o',facecolors='none',edgecolors='black',linewidths=.7)
        ax.set_xticks(range(4),['MHC-I','MHC-II','Cross-pres.','IFN-gamma'],rotation=45,ha='right')
        ax.set_yticks(range(len(labels)),labels);ax.tick_params(axis='y',labelsize=6.6)
        panel(ax,'G','Candidate genes and antigen-presentation pathways')
        cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.04);cb.set_label('-log10 module q (q < 0.05)',fontsize=6.5)
        ax.text(.0,-.29,'Circle: member in Reactome 2022',transform=ax.transAxes,fontsize=6.4)
        # H-K. Group pathway curves.
        pp=base/'disease_pathways/pathway_pseudotime';curves=pd.read_csv(pp/'antigen_pathway_fitted_curves.csv')
        primary=pd.read_csv(pp/'antigen_pathway_disease_pseudotime_tests.csv').set_index('pathway')
        sens=pd.read_csv(pp/'antigen_source_curve_sensitivity.csv').set_index('pathway')
        for i,(pid,name) in enumerate([('R-HSA-983170','MHC-I peptide loading'),('R-HSA-2132295','MHC-II presentation'),
         ('R-HSA-1236975','Cross-presentation'),('R-HSA-877300','IFN-gamma signaling')]):
         ax=fig.add_subplot(gs[3,3*i:3*(i+1)]);path=next(p for p in primary.index if pid in p)
         for group,color in colors.items():
          d=curves[curves.pathway.eq(path)&curves.status.eq(group)]
          ax.plot(d.pseudotime,d.fitted_score,color=color,label=group,lw=.9)
          ax.fill_between(d.pseudotime,d.ci_low,d.ci_high,color=color,alpha=.10,lw=0)
         ax.set(xlabel='Pseudotime',ylabel='Relative pathway score');panel(ax,chr(ord('H')+i),name)
         ax.text(.02,.98,'Shape q='+format(primary.loc[path,'FDR_shape'],'.3g')+'\nSource-curve q='+format(sens.loc[path,'FDR_shape'],'.3g'),transform=ax.transAxes,va='top',fontsize=6.4)
         if i==0:ax.legend(frameon=False,fontsize=6.3,loc='lower left')
        for a in fig.axes:
         a.set_rasterization_zorder(None)
         for c in a.collections:c.set_rasterized(False)
        svg=out/'Myeloid_COPD_Composite.svg'
        fig.savefig(svg,format='svg',metadata={'Title':'Myeloid COPD composite','Description':'Panels A-K, with original scPagwas PDF, GSEA and candidate-pathway associations.'})
        # Import the specified PDF through its vector SVG representation; preserve all source paths.
        S='http://www.w3.org/2000/svg';X='http://www.w3.org/1999/xlink'
        ET.register_namespace('',S);ET.register_namespace('xlink',X)
        tree=ET.parse(svg);destination=tree.getroot();source=ET.parse(pdf_svg).getroot()
        source_ids={e.attrib['id'] for e in source.iter() if 'id' in e.attrib}
        idmap={v:'panelA_pdf_'+v for v in source_ids}
        for elem in source.iter():
         for k,v in list(elem.attrib.items()):
          if k=='id':elem.set(k,idmap[v]);continue
          if v.startswith('#') and v[1:] in idmap:elem.set(k,'#'+idmap[v[1:]])
          else:
           for old,new in idmap.items():v=v.replace('url(#'+old+')','url(#'+new+')')
           elem.set(k,v)
        pos=axA.get_position();W,H=fig.get_size_inches()*72
        left=pos.x0*W;top=(1-pos.y1)*H
        width=pos.width*W;height=pos.height*H
        scale=min(width/504,height/504);left+=(width-504*scale)/2;top+=(height-504*scale)/2
        wrapper=ET.Element('{'+S+'}g',{'id':'panelA_original_PDF','transform':f'translate({left} {top}) scale({scale})'})
        for element in list(source):wrapper.append(copy.deepcopy(element))
        destination.append(wrapper);tree.write(svg,encoding='utf-8',xml_declaration=True)
        plt.close(fig)

if __name__ == "__main__":
    main()
