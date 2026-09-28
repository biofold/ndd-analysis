import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt

def compute_roc_pr(filename, predictor_col, outcome_col, output_file=None, pos_label=1,
                    save_pdf=True, save_png=False, save_stats=True):
    """
    Compute the ROC curve/AUC and Precision-Recall curve/AUC (average
    precision) for a continuous or ordinal predictor against a binary (0/1)
    outcome, and plot both curves side by side.

    Mirrors the file-reading, column-indexing, and plotting conventions of
    violin.py / bar.py / ks-test.py / fisher-test.py (1-based column
    indices, white background, PDF-friendly fonts, full plot borders), so
    it can be dropped into the same MOE-score validation pipeline.

    Args:
        filename: Path to input file (CSV, Excel, or TSV)
        predictor_col: 1-based column index for the predictor score (e.g. MOE_score)
        outcome_col: 1-based column index for the binary (0/1) outcome
        output_file: Output file name without extension (optional)
        pos_label: Which value of the outcome column counts as the positive class (default: 1)
        save_pdf: Whether to save as PDF (default: True)
        save_png: Whether to save as PNG (default: False)
        save_stats: Whether to save the curve points and summary AUCs to file (default: True)
    """
    try:
        from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
    except ImportError:
        print("Error: scikit-learn is required (pip install scikit-learn).")
        return None

    file_ext = os.path.splitext(filename)[1].lower()

    try:
        if file_ext == '.csv':
            df = pd.read_csv(filename)
        elif file_ext in ['.xlsx', '.xls']:
            df = pd.read_excel(filename)
        elif file_ext in ['.tsv', '.txt']:
            df = pd.read_csv(filename, sep='\t')
        else:
            df = pd.read_csv(filename, sep=None, engine='python')
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

    pred_idx = predictor_col - 1
    outcome_idx = outcome_col - 1

    if pred_idx >= len(df.columns) or outcome_idx >= len(df.columns):
        print(f"Error: Column index out of range. File has {len(df.columns)} columns.")
        return None

    pred_col_name = df.columns[pred_idx]
    outcome_col_name = df.columns[outcome_idx]

    print(f"Predictor column: '{pred_col_name}' (column {predictor_col})")
    print(f"Outcome column (binary): '{outcome_col_name}' (column {outcome_col})")

    df_clean = df[[pred_col_name, outcome_col_name]].dropna()
    df_clean[pred_col_name] = pd.to_numeric(df_clean[pred_col_name], errors='coerce')
    df_clean[outcome_col_name] = pd.to_numeric(df_clean[outcome_col_name], errors='coerce')
    df_clean = df_clean.dropna()

    y = (df_clean[outcome_col_name].values == pos_label).astype(int)
    score = df_clean[pred_col_name].values

    n = len(y)
    n_pos = int(y.sum())
    prevalence = n_pos / n if n else float('nan')
    print(f"\nTotal data points: {n}")
    print(f"Positives: {n_pos} ({prevalence:.1%})")

    fpr, tpr, _ = roc_curve(y, score)
    roc_auc = auc(fpr, tpr)
    precision, recall, _ = precision_recall_curve(y, score)
    pr_auc = average_precision_score(y, score)

    print(f"\nROC-AUC: {roc_auc:.4f}")
    print(f"PR-AUC (average precision): {pr_auc:.4f}")
    print(f"Baseline (prevalence): {prevalence:.4f}")

    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), facecolor='white')

    ax0 = axes[0]
    ax0.set_facecolor('white')
    ax0.plot(fpr, tpr, color='#3B6FA0', linewidth=1.8, label=f'AUC = {roc_auc:.3f}')
    ax0.plot([0, 1], [0, 1], color='#999999', linestyle='--', linewidth=1)
    ax0.set_xlabel('False Positive Rate', fontsize=13, labelpad=10)
    ax0.set_ylabel('True Positive Rate', fontsize=13, labelpad=10)
    ax0.set_title(f'ROC — {outcome_col_name} vs {pred_col_name}', fontsize=13, loc='left')
    ax0.legend(loc='lower right', frameon=False, fontsize=9)
    for spine in ax0.spines.values():
        spine.set_visible(True)
        spine.set_color('#333333')
        spine.set_linewidth(0.8)

    ax1 = axes[1]
    ax1.set_facecolor('white')
    ax1.plot(recall, precision, color='#C0504D', linewidth=1.8, label=f'AUC = {pr_auc:.3f}')
    ax1.axhline(prevalence, color='#999999', linestyle='--', linewidth=1, label=f'baseline = {prevalence:.3f}')
    ax1.set_xlabel('Recall', fontsize=13, labelpad=10)
    ax1.set_ylabel('Precision', fontsize=13, labelpad=10)
    ax1.set_title(f'PR — {outcome_col_name} vs {pred_col_name}', fontsize=13, loc='left')
    ax1.legend(loc='upper right', frameon=False, fontsize=9)
    for spine in ax1.spines.values():
        spine.set_visible(True)
        spine.set_color('#333333')
        spine.set_linewidth(0.8)

    plt.tight_layout()

    if output_file is None:
        base_filename = f"roc_pr_{pred_col_name}_vs_{outcome_col_name}"
    else:
        base_filename = os.path.splitext(output_file)[0]

    saved_files = []

    if save_pdf:
        pdf_file = f"{base_filename}.pdf"
        plt.savefig(pdf_file, format='pdf', bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"\n✓ PDF saved as: {pdf_file}")
        saved_files.append(pdf_file)

    if save_png:
        png_file = f"{base_filename}.png"
        plt.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"✓ PNG saved as: {png_file}")
        saved_files.append(png_file)

    summary_df = pd.DataFrame([{
        'Predictor': pred_col_name,
        'Outcome': outcome_col_name,
        'N': n,
        'N_positive': n_pos,
        'Prevalence': prevalence,
        'ROC_AUC': roc_auc,
        'PR_AUC': pr_auc,
    }])

    if save_stats:
        stats_file = f"{base_filename}_summary.tsv"
        summary_df.to_csv(stats_file, sep='\t', index=False)
        print(f"✓ Summary AUC table saved as: {stats_file}")
        saved_files.append(stats_file)

        curve_file = f"{base_filename}_curves.tsv"
        roc_df = pd.DataFrame({'curve': 'ROC', 'x_fpr_or_recall': fpr, 'y_tpr_or_precision': tpr})
        pr_df = pd.DataFrame({'curve': 'PR', 'x_fpr_or_recall': recall, 'y_tpr_or_precision': precision})
        pd.concat([roc_df, pr_df], ignore_index=True).to_csv(curve_file, sep='\t', index=False)
        print(f"✓ Curve points saved as: {curve_file}")
        saved_files.append(curve_file)

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(summary_df.to_string(index=False))

    return summary_df

def main():
    """Main function to handle command line arguments"""

    # Optional --png flag: also write a PNG alongside the PDF. The manuscript
    # figures are vector PDFs, but the response-to-reviewers document embeds
    # raster copies, so both are needed from the same invocation.
    save_png = "--png" in sys.argv
    if save_png:
        sys.argv = [a for a in sys.argv if a != "--png"]


    if len(sys.argv) >= 4:
        filename = sys.argv[1]
        predictor_col = int(sys.argv[2])
        outcome_col = int(sys.argv[3])
        output_file = sys.argv[4] if len(sys.argv) > 4 else None
        pos_label = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    else:
        print("ROC / Precision-Recall AUC Calculator")
        print("-" * 60)

        filename = input("Enter the input file path: ").strip()

        if not os.path.exists(filename):
            print(f"Error: File '{filename}' does not exist.")
            return

        try:
            predictor_col = int(input("Enter the predictor column number (1-based): "))
            outcome_col = int(input("Enter the binary outcome column number (1-based): "))
        except ValueError:
            print("Error: Column numbers must be integers.")
            return

        output_file = input("Enter output file name (or press Enter for default): ").strip()
        if not output_file:
            output_file = None
        pos_label = 1

    compute_roc_pr(filename, predictor_col, outcome_col, output_file, pos_label,
                   save_pdf=True, save_png=save_png, save_stats=True)

if __name__ == "__main__":
    required_packages = ['pandas', 'numpy', 'matplotlib', 'scikit-learn']

    import subprocess

    for package in required_packages:
        try:
            __import__(package if package != 'scikit-learn' else 'sklearn')
        except ImportError:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])

    main()
