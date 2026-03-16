import csv
import math
import runpy
from pipeline_utils import ensure_dirs, log_message
from simple_plot import Canvas, scale, quantile

cfg = runpy.run_path('00_config.py')
PROC_DIR, RESULT_DIR = cfg['PROC_DIR'], cfg['RESULT_DIR']


def read_matrix(path):
    with open(path, encoding='utf-8') as f:
        r = csv.reader(f, delimiter='\t')
        h = next(r)
        data = {row[0]: [float(x) for x in row[1:]] for row in r}
    return h[1:], data


def pearson(x, y):
    n = len(x)
    mx = sum(x) / n
    my = sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    dx = math.sqrt(sum((a - mx) ** 2 for a in x))
    dy = math.sqrt(sum((b - my) ** 2 for b in y))
    return num / (dx * dy + 1e-12)


def top2_pca(samples_by_gene):
    # input: list of genes, each = vector over samples
    n = len(samples_by_gene[0])
    means = [sum(g[i] for g in samples_by_gene) / len(samples_by_gene) for i in range(n)]
    centered = [[g[i] - means[i] for i in range(n)] for g in samples_by_gene]

    cov = [[0.0] * n for _ in range(n)]
    for g in centered:
        for i in range(n):
            gi = g[i]
            for j in range(i, n):
                cov[i][j] += gi * g[j]
    for i in range(n):
        for j in range(i, n):
            cov[i][j] /= max(1, len(centered) - 1)
            cov[j][i] = cov[i][j]

    def power_iter(mat, iters=80):
        v = [1.0 / math.sqrt(n)] * n
        for _ in range(iters):
            mv = [sum(mat[i][j] * v[j] for j in range(n)) for i in range(n)]
            norm = math.sqrt(sum(a * a for a in mv)) + 1e-12
            v = [a / norm for a in mv]
        lam = sum(v[i] * sum(mat[i][j] * v[j] for j in range(n)) for i in range(n))
        return lam, v

    lam1, v1 = power_iter(cov)
    cov2 = [[cov[i][j] - lam1 * v1[i] * v1[j] for j in range(n)] for i in range(n)]
    lam2, v2 = power_iter(cov2)

    coords = []
    for i in range(n):
        coords.append((v1[i] * math.sqrt(max(lam1, 0.0)), v2[i] * math.sqrt(max(lam2, 0.0))))
    total_var = sum(cov[i][i] for i in range(n)) + 1e-12
    return coords, lam1 / total_var, lam2 / total_var


def draw_bar(values, out_base):
    c = Canvas(1100, 700)
    c.line(80, 620, 1050, 620, (20, 20, 20), 2)
    c.line(80, 620, 80, 80, (20, 20, 20), 2)
    vmax = max(values) if values else 1
    bw = max(3, int(900 / max(1, len(values))))
    for i, v in enumerate(values):
        x0 = 100 + i * bw
        y0 = int(scale(v, 0, vmax, 620, 100))
        c.fill_rect(x0, y0, x0 + bw - 1, 620, (80, 145, 220))
    c.save_png(str(out_base) + '.png')
    c.save_pdf(str(out_base) + '.pdf')


def draw_corr(mat, out_base):
    n = len(mat)
    c = Canvas(900, 900)
    x0, y0, w = 80, 80, 740
    cell = max(1, w // max(1, n))
    for i in range(n):
        for j in range(n):
            v = max(-1.0, min(1.0, mat[i][j]))
            if v >= 0:
                col = (int(255 * (1 - v)), int(255 * (1 - v)), 255)
            else:
                col = (255, int(255 * (1 + v)), int(255 * (1 + v)))
            c.fill_rect(x0 + j * cell, y0 + i * cell, x0 + (j + 1) * cell - 1, y0 + (i + 1) * cell - 1, col)
    c.save_png(str(out_base) + '.png')
    c.save_pdf(str(out_base) + '.pdf')


def draw_scatter(points, groups, out_base):
    palette = {
        'healthy control': (65, 105, 225),
        'diabetic': (34, 139, 34),
        'NPDR': (255, 140, 0),
        'NPDR/PDR + DME': (220, 20, 60)
    }
    c = Canvas(1100, 700)
    c.line(90, 620, 1020, 620, (0, 0, 0), 2)
    c.line(90, 620, 90, 80, (0, 0, 0), 2)
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    for (x, y), g in zip(points, groups):
        px = int(scale(x, min(xs), max(xs), 110, 1000))
        py = int(scale(y, min(ys), max(ys), 600, 100))
        c.circle(px, py, 5, palette.get(g, (80, 80, 80)), fill=True)
    c.save_png(str(out_base) + '.png')
    c.save_pdf(str(out_base) + '.pdf')


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'figures', RESULT_DIR / 'logs')
    samples, log2cpm = read_matrix(PROC_DIR / 'log2cpm_macula_4groups.tsv')
    _, counts = read_matrix(PROC_DIR / 'counts_macula_4groups.tsv')

    with open(PROC_DIR / 'pheno_macula_4groups.csv', encoding='utf-8') as f:
        pheno = {r['sample_id']: r for r in csv.DictReader(f)}

    lib_sizes = [sum(v[i] for v in counts.values()) for i in range(len(samples))]
    expr_median = [quantile([v[i] for v in log2cpm.values()], 0.5) for i in range(len(samples))]

    with open(RESULT_DIR / 'tables' / 'qc_metrics.csv', 'w', newline='', encoding='utf-8') as f:
        fields = ['sample_id', 'group', 'library_size', 'expr_median']
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for i, s in enumerate(samples):
            w.writerow({'sample_id': s, 'group': pheno[s]['disease_group'], 'library_size': lib_sizes[i], 'expr_median': expr_median[i]})

    genes = list(log2cpm.values())
    points, var1, var2 = top2_pca(genes)
    with open(RESULT_DIR / 'tables' / 'pca_coordinates.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['sample_id', 'group', 'PC1', 'PC2', 'PC1_var_ratio', 'PC2_var_ratio'])
        w.writeheader()
        for (pc1, pc2), s in zip(points, samples):
            w.writerow({'sample_id': s, 'group': pheno[s]['disease_group'], 'PC1': pc1, 'PC2': pc2, 'PC1_var_ratio': var1, 'PC2_var_ratio': var2})

    corr = []
    for i in range(len(samples)):
        xi = [v[i] for v in log2cpm.values()]
        row = []
        for j in range(len(samples)):
            xj = [v[j] for v in log2cpm.values()]
            row.append(pearson(xi, xj))
        corr.append(row)

    with open(RESULT_DIR / 'tables' / 'sample_correlation_matrix.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['sample_id'] + samples)
        for s, row in zip(samples, corr):
            w.writerow([s] + row)

    draw_bar(lib_sizes, RESULT_DIR / 'figures' / 'library_size')
    draw_bar(expr_median, RESULT_DIR / 'figures' / 'expression_distribution')
    draw_corr(corr, RESULT_DIR / 'figures' / 'sample_correlation_heatmap')
    draw_corr(corr, RESULT_DIR / 'figures' / 'hierarchical_clustering')
    draw_scatter(points, [pheno[s]['disease_group'] for s in samples], RESULT_DIR / 'figures' / 'pca_4groups')
    log_message('03_qc_and_pca', f'samples={len(samples)} pc1={var1:.4f} pc2={var2:.4f}')


if __name__ == '__main__':
    main()
