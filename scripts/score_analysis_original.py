import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# 设置中文显示
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def load_and_explore_data(filepath):
    """加载数据并进行初步探索"""
    print("=" * 50)
    print("Loading and Exploring Data")
    print("=" * 50)
    
    # 读取CSV文件
    df = pd.read_csv(filepath)
    
    print(f"\nDataset Shape: {df.shape}")
    print(f"Total Records: {len(df):,}")
    print(f"Number of Columns: {df.shape[1]}")
    
    # 基本信息
    print("\nColumn Information:")
    print(df.info())
    
    # 前几行数据
    print("\nFirst 5 rows:")
    print(df.head())
    
    # 数据类型
    print("\nData Types:")
    print(df.dtypes)
    
    # 缺失值检查
    print("\nMissing Values:")
    missing = df.isnull().sum()
    if missing.any():
        print(missing[missing > 0])
    else:
        print("No missing values found")
    
    return df

def analyze_score_distribution(df):
    """分析结合能评分分布"""
    print("\n" + "=" * 50)
    print("Binding Score Distribution Analysis")
    print("=" * 50)
    
    # KPNA4评分统计
    print("\nKPNA4 Score Statistics:")
    print(f"  Mean: {df['Score_KPNA4'].mean():.3f}")
    print(f"  Median: {df['Score_KPNA4'].median():.3f}")
    print(f"  Std Dev: {df['Score_KPNA4'].std():.3f}")
    print(f"  Min: {df['Score_KPNA4'].min():.3f}")
    print(f"  Max: {df['Score_KPNA4'].max():.3f}")
    print(f"  25th Percentile: {df['Score_KPNA4'].quantile(0.25):.3f}")
    print(f"  75th Percentile: {df['Score_KPNA4'].quantile(0.75):.3f}")
    
    # MC1R评分统计
    print("\nMC1R Score Statistics:")
    print(f"  Mean: {df['Score_MC1R'].mean():.3f}")
    print(f"  Median: {df['Score_MC1R'].median():.3f}")
    print(f"  Std Dev: {df['Score_MC1R'].std():.3f}")
    print(f"  Min: {df['Score_MC1R'].min():.3f}")
    print(f"  Max: {df['Score_MC1R'].max():.3f}")
    print(f"  25th Percentile: {df['Score_MC1R'].quantile(0.25):.3f}")
    print(f"  75th Percentile: {df['Score_MC1R'].quantile(0.75):.3f}")
    
    # 相关性分析
    correlation = df['Score_KPNA4'].corr(df['Score_MC1R'])
    print(f"\nCorrelation between KPNA4 and MC1R scores: {correlation:.3f}")
    
    # 创建可视化
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # KPNA4分布直方图
    axes[0, 0].hist(df['Score_KPNA4'], bins=50, edgecolor='black', alpha=0.7)
    axes[0, 0].set_title('KPNA4 Score Distribution')
    axes[0, 0].set_xlabel('Score')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].axvline(df['Score_KPNA4'].mean(), color='red', linestyle='--', label='Mean')
    axes[0, 0].axvline(df['Score_KPNA4'].median(), color='green', linestyle='--', label='Median')
    axes[0, 0].legend()
    
    # MC1R分布直方图
    axes[0, 1].hist(df['Score_MC1R'], bins=50, edgecolor='black', alpha=0.7, color='orange')
    axes[0, 1].set_title('MC1R Score Distribution')
    axes[0, 1].set_xlabel('Score')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].axvline(df['Score_MC1R'].mean(), color='red', linestyle='--', label='Mean')
    axes[0, 1].axvline(df['Score_MC1R'].median(), color='green', linestyle='--', label='Median')
    axes[0, 1].legend()
    
    # 散点图
    axes[0, 2].scatter(df['Score_KPNA4'], df['Score_MC1R'], alpha=0.5, s=1)
    axes[0, 2].set_title(f'KPNA4 vs MC1R Scores (r={correlation:.3f})')
    axes[0, 2].set_xlabel('KPNA4 Score')
    axes[0, 2].set_ylabel('MC1R Score')
    
    # 箱线图比较
    box_data = pd.DataFrame({
        'KPNA4': df['Score_KPNA4'],
        'MC1R': df['Score_MC1R']
    })
    box_data_melted = box_data.melt(var_name='Target', value_name='Score')
    axes[1, 0].boxplot([df['Score_KPNA4'], df['Score_MC1R']], labels=['KPNA4', 'MC1R'])
    axes[1, 0].set_title('Score Distribution Comparison')
    axes[1, 0].set_ylabel('Score')
    
    # QQ图 - KPNA4
    stats.probplot(df['Score_KPNA4'], dist="norm", plot=axes[1, 1])
    axes[1, 1].set_title('Q-Q Plot - KPNA4')
    
    # QQ图 - MC1R
    stats.probplot(df['Score_MC1R'], dist="norm", plot=axes[1, 2])
    axes[1, 2].set_title('Q-Q Plot - MC1R')
    
    plt.tight_layout()
    plt.savefig('score_distribution_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def analyze_peptide_structure(df):
    """分析多肽结构特征"""
    print("\n" + "=" * 50)
    print("Peptide Structure Analysis")
    print("=" * 50)
    
    # Root分析
    print(f"\nUnique Roots: {df['Root'].nunique():,}")
    root_counts = df['Root'].value_counts()
    print(f"Most Common Roots (Top 10):")
    for i, (root, count) in enumerate(root_counts.head(10).items(), 1):
        print(f"  {i}. {root}: {count:,} peptides")
    
    # Part分析
    print(f"\nUnique Parts: {df['Part'].nunique():,}")
    part_counts = df['Part'].value_counts()
    print(f"Most Common Parts (Top 10):")
    for i, (part, count) in enumerate(part_counts.head(10).items(), 1):
        print(f"  {i}. {part}: {count:,} peptides")
    
    # 配体序列长度分析
    df['Ligand_Length'] = df['Ligand'].str.len()
    print(f"\nLigand Sequence Length Statistics:")
    print(f"  Mean Length: {df['Ligand_Length'].mean():.1f}")
    print(f"  Median Length: {df['Ligand_Length'].median():.0f}")
    print(f"  Min Length: {df['Ligand_Length'].min()}")
    print(f"  Max Length: {df['Ligand_Length'].max()}")
    
    # 氨基酸组成分析
    print("\nAmino Acid Composition Analysis:")
    aa_counts = {}
    for seq in df['Ligand']:
        for aa in seq:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
    
    aa_df = pd.DataFrame(list(aa_counts.items()), columns=['AA', 'Count'])
    aa_df = aa_df.sort_values('Count', ascending=False)
    
    print("Top 10 Most Frequent Amino Acids:")
    for i, row in aa_df.head(10).iterrows():
        percentage = (row['Count'] / sum(aa_counts.values())) * 100
        print(f"  {row['AA']}: {row['Count']:,} ({percentage:.2f}%)")

def identify_top_binders(df, top_n=100):
    """识别最优结合多肽"""
    print("\n" + "=" * 50)
    print(f"Top {top_n} Binding Peptides Analysis")
    print("=" * 50)
    
    # KPNA4最优结合肽
    top_kpna4 = df.nsmallest(top_n, 'Score_KPNA4')
    print(f"\nTop {top_n} KPNA4 Binders:")
    print(f"  Score Range: [{top_kpna4['Score_KPNA4'].max():.3f}, {top_kpna4['Score_KPNA4'].min():.3f}]")
    print(f"  Mean Score: {top_kpna4['Score_KPNA4'].mean():.3f}")
    
    # MC1R最优结合肽
    top_mc1r = df.nsmallest(top_n, 'Score_MC1R')
    print(f"\nTop {top_n} MC1R Binders:")
    print(f"  Score Range: [{top_mc1r['Score_MC1R'].max():.3f}, {top_mc1r['Score_MC1R'].min():.3f}]")
    print(f"  Mean Score: {top_mc1r['Score_MC1R'].mean():.3f}")
    
    # 双靶标优秀结合肽
    df['Combined_Score'] = df['Score_KPNA4'] + df['Score_MC1R']
    top_both = df.nsmallest(top_n, 'Combined_Score')
    print(f"\nTop {top_n} Dual-Target Binders (Combined Score):")
    print(f"  Combined Score Range: [{top_both['Combined_Score'].max():.3f}, {top_both['Combined_Score'].min():.3f}]")
    
    # 选择性分析
    df['Selectivity_KPNA4'] = df['Score_MC1R'] - df['Score_KPNA4']  # 正值表示对KPNA4选择性更高
    df['Selectivity_MC1R'] = df['Score_KPNA4'] - df['Score_MC1R']   # 正值表示对MC1R选择性更高
    
    selective_kpna4 = df.nlargest(top_n, 'Selectivity_KPNA4')
    selective_mc1r = df.nlargest(top_n, 'Selectivity_MC1R')
    
    print(f"\nTop {top_n} KPNA4-Selective Peptides:")
    print(f"  Mean Selectivity: {selective_kpna4['Selectivity_KPNA4'].mean():.3f}")
    
    print(f"\nTop {top_n} MC1R-Selective Peptides:")
    print(f"  Mean Selectivity: {selective_mc1r['Selectivity_MC1R'].mean():.3f}")
    
    return top_kpna4, top_mc1r, top_both

def analyze_sequence_patterns(df, top_binders):
    """分析优秀结合肽的序列模式"""
    print("\n" + "=" * 50)
    print("Sequence Pattern Analysis in Top Binders")
    print("=" * 50)
    
    # 分析顶部结合肽的Root和Part模式
    print("\nMost Common Roots in Top Binders:")
    root_pattern = top_binders['Root'].value_counts().head(10)
    for i, (root, count) in enumerate(root_pattern.items(), 1):
        percentage = (count / len(top_binders)) * 100
        print(f"  {i}. {root}: {count} ({percentage:.1f}%)")
    
    print("\nMost Common Parts in Top Binders:")
    part_pattern = top_binders['Part'].value_counts().head(10)
    for i, (part, count) in enumerate(part_pattern.items(), 1):
        percentage = (count / len(top_binders)) * 100
        print(f"  {i}. {part}: {count} ({percentage:.1f}%)")
    
    # 分析序列长度偏好
    top_length_mean = top_binders['Ligand'].str.len().mean()
    all_length_mean = df['Ligand'].str.len().mean()
    print(f"\nSequence Length Preference:")
    print(f"  Top Binders Mean Length: {top_length_mean:.1f}")
    print(f"  All Peptides Mean Length: {all_length_mean:.1f}")
    print(f"  Difference: {top_length_mean - all_length_mean:+.1f}")

def create_heatmap_analysis(df, sample_size=1000):
    """创建热图分析"""
    print("\n" + "=" * 50)
    print("Creating Heatmap Analysis")
    print("=" * 50)
    
    # 采样数据用于热图（如果数据太大）
    if len(df) > sample_size:
        df_sample = df.sample(n=sample_size, random_state=42)
        print(f"Sampling {sample_size} peptides for visualization")
    else:
        df_sample = df
    
    # 按Root和Part分组的平均得分
    pivot_kpna4 = df.pivot_table(values='Score_KPNA4', 
                                  index=df['Root'].value_counts().head(20).index,
                                  columns=df['Part'].value_counts().head(20).index,
                                  aggfunc='mean')
    
    pivot_mc1r = df.pivot_table(values='Score_MC1R',
                                index=df['Root'].value_counts().head(20).index,
                                columns=df['Part'].value_counts().head(20).index,
                                aggfunc='mean')
    
    # 创建热图
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    
    sns.heatmap(pivot_kpna4, cmap='coolwarm', center=0, 
                ax=axes[0], cbar_kws={'label': 'Mean KPNA4 Score'})
    axes[0].set_title('KPNA4 Binding Score Heatmap\n(Root vs Part)')
    axes[0].set_xlabel('Part')
    axes[0].set_ylabel('Root')
    
    sns.heatmap(pivot_mc1r, cmap='coolwarm', center=0,
                ax=axes[1], cbar_kws={'label': 'Mean MC1R Score'})
    axes[1].set_title('MC1R Binding Score Heatmap\n(Root vs Part)')
    axes[1].set_xlabel('Part')
    axes[1].set_ylabel('Root')
    
    plt.tight_layout()
    plt.savefig('binding_score_heatmap.png', dpi=300, bbox_inches='tight')
    plt.show()

def save_top_peptides(top_kpna4, top_mc1r, top_both, output_prefix='top_peptides'):
    """保存最优结合肽到文件"""
    print("\n" + "=" * 50)
    print("Saving Top Peptides to Files")
    print("=" * 50)
    
    # 保存KPNA4最优结合肽
    top_kpna4_file = f"{output_prefix}_kpna4.csv"
    top_kpna4.to_csv(top_kpna4_file, index=False)
    print(f"Saved top KPNA4 binders to: {top_kpna4_file}")
    
    # 保存MC1R最优结合肽
    top_mc1r_file = f"{output_prefix}_mc1r.csv"
    top_mc1r.to_csv(top_mc1r_file, index=False)
    print(f"Saved top MC1R binders to: {top_mc1r_file}")
    
    # 保存双靶标最优结合肽
    top_both_file = f"{output_prefix}_dual_target.csv"
    top_both.to_csv(top_both_file, index=False)
    print(f"Saved top dual-target binders to: {top_both_file}")

def main(filepath='merged_scores_sp_id.csv'):
    """主函数：执行完整的多肽结合能分析"""
    print("\n" + "="*60)
    print(" PEPTIDE BINDING ENERGY LANDSCAPE ANALYSIS ")
    print("="*60)
    
    try:
        # 1. 加载和探索数据
        df = load_and_explore_data(filepath)
        
        # 2. 分析得分分布
        analyze_score_distribution(df)
        
        # 3. 分析多肽结构
        analyze_peptide_structure(df)
        
        # 4. 识别最优结合肽
        top_kpna4, top_mc1r, top_both = identify_top_binders(df, top_n=100)
        
        # 5. 分析序列模式
        analyze_sequence_patterns(df, top_both)
        
        # 6. 创建热图分析
        create_heatmap_analysis(df)
        
        # 7. 保存结果
        save_top_peptides(top_kpna4, top_mc1r, top_both)
        
        # 8. 生成总结报告
        print("\n" + "="*60)
        print(" ANALYSIS SUMMARY ")
        print("="*60)
        
        print(f"\nDataset Overview:")
        print(f"  Total Peptides Analyzed: {len(df):,}")
        print(f"  Unique Roots: {df['Root'].nunique():,}")
        print(f"  Unique Parts: {df['Part'].nunique():,}")
        
        print(f"\nBinding Score Ranges:")
        print(f"  KPNA4: [{df['Score_KPNA4'].min():.3f}, {df['Score_KPNA4'].max():.3f}]")
        print(f"  MC1R: [{df['Score_MC1R'].min():.3f}, {df['Score_MC1R'].max():.3f}]")
        
        print(f"\nTop Performers:")
        print(f"  Best KPNA4 Binder: {top_kpna4.iloc[0]['Ligand']} (Score: {top_kpna4.iloc[0]['Score_KPNA4']:.3f})")
        print(f"  Best MC1R Binder: {top_mc1r.iloc[0]['Ligand']} (Score: {top_mc1r.iloc[0]['Score_MC1R']:.3f})")
        print(f"  Best Dual-Target Binder: {top_both.iloc[0]['Ligand']}")
        print(f"    KPNA4 Score: {top_both.iloc[0]['Score_KPNA4']:.3f}")
        print(f"    MC1R Score: {top_both.iloc[0]['Score_MC1R']:.3f}")
        
        print("\n" + "="*60)
        print(" Analysis Complete! ")
        print("="*60)
        
    except Exception as e:
        print(f"\nError during analysis: {str(e)}")
        print("Please check the file path and data format.")
        raise

if __name__ == "__main__":
    # 运行主分析
    main('merged_scores_sp_id.csv')
