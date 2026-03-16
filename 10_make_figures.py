import runpy
from pipeline_utils import ensure_dirs, log_message
cfg=runpy.run_path('00_config.py')
RESULT_DIR=cfg['RESULT_DIR']

def touch_fig(name):
    for ext in cfg['FIG_FORMATS']:
        (RESULT_DIR/'figures'/f'{name}.{ext}').write_text(f'placeholder {name} dpi={cfg["FIG_DPI"]}',encoding='utf-8')

def main():
    ensure_dirs(RESULT_DIR/'figures', RESULT_DIR/'logs')
    figs=['Figure1_workflow','Figure1_sample_composition','Figure1_pca_macula_4groups','Figure2_inflammation_ssgsea_4groups','Figure2_primary_volcano','Figure2_deg_inflammatory_heatmap','Figure2_overlap_primarydeg_inflammation','Figure3_lasso_path','Figure3_cv_curve','Figure3_selected_gene_expression_4groups','Figure3_individual_roc','Figure3_combined_roc','Figure4_selected_gene_gsea_heatmap','Figure4_recurrent_pathway_dotplot','Figure4_gene_pathway_network_optional','Figure5_immune_score_comparison','Figure5_gene_immune_corr_heatmap','Supp_QC','Supp_clustering','Supp_expanded_deg_heatmap','Supp_all_inflammatory_core_heatmap','Supp_single_gene_gsea','Supp_gene_immune_scatter']
    for f in figs: touch_fig(f)
    log_message('10_make_figures',f'generated={len(figs)} placeholders')
if __name__=='__main__': main()
