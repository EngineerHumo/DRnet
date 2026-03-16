import base64
import runpy
from pipeline_utils import ensure_dirs, log_message

cfg = runpy.run_path('00_config.py')
RESULT_DIR = cfg['RESULT_DIR']

# 1x1 transparent PNG
PNG_BYTES = base64.b64decode(
    'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO7Z0r8AAAAASUVORK5CYII='
)
PDF_BYTES = b'%PDF-1.1\n1 0 obj<< /Type /Catalog /Pages 2 0 R >>endobj\n2 0 obj<< /Type /Pages /Kids [3 0 R] /Count 1 >>endobj\n3 0 obj<< /Type /Page /Parent 2 0 R /MediaBox [0 0 100 100] /Contents 4 0 R >>endobj\n4 0 obj<< /Length 35 >>stream\nBT /F1 12 Tf 10 50 Td (placeholder) Tj ET\nendstream\nendobj\nxref\n0 5\n0000000000 65535 f \n0000000010 00000 n \n0000000061 00000 n \n0000000118 00000 n \n0000000210 00000 n \ntrailer<< /Root 1 0 R /Size 5 >>\nstartxref\n295\n%%EOF\n'


def touch_fig(name):
    png = RESULT_DIR / 'figures' / f'{name}.png'
    pdf = RESULT_DIR / 'figures' / f'{name}.pdf'
    png.write_bytes(PNG_BYTES)
    pdf.write_bytes(PDF_BYTES)


def main():
    ensure_dirs(RESULT_DIR / 'figures', RESULT_DIR / 'logs')
    figs = [
        'Figure1_workflow', 'Figure1_sample_composition', 'Figure1_pca_macula_4groups',
        'Figure2_inflammation_ssgsea_4groups', 'Figure2_primary_volcano', 'Figure2_deg_inflammatory_heatmap',
        'Figure2_overlap_primarydeg_inflammation', 'Figure3_lasso_path', 'Figure3_cv_curve',
        'Figure3_selected_gene_expression_4groups', 'Figure3_individual_roc', 'Figure3_combined_roc',
        'Figure4_selected_gene_gsea_heatmap', 'Figure4_recurrent_pathway_dotplot', 'Figure4_gene_pathway_network_optional',
        'Figure5_immune_score_comparison', 'Figure5_gene_immune_corr_heatmap', 'Supp_QC', 'Supp_clustering',
        'Supp_expanded_deg_heatmap', 'Supp_all_inflammatory_core_heatmap', 'Supp_single_gene_gsea', 'Supp_gene_immune_scatter'
    ]
    for f in figs:
        touch_fig(f)
    log_message('10_make_figures', f'generated={len(figs)} loadable placeholder images (PNG/PDF)')


if __name__ == '__main__':
    main()
