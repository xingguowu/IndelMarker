# InDelMarker：策略 1/2/3 统一 Snakemake 工作流

先用随工作流提供的 Conda 环境文件安装全部依赖。它将 Python 固定在 3.10，因为 SyRI 1.8.2 需要该版本；Snakemake 保持在 8.x 以与其共存。Bioconda 当前同时提供 SyRI 1.8.2 的 `osx-arm64` 与 `linux-aarch64` 构建，[SyRI recipe](https://bioconda.github.io/recipes/syri/README.html)；Primer3 和 Minimap2 也有 `osx-arm64` 构建，[Primer3 recipe](https://bioconda.github.io/recipes/primer3/README.html)、[Minimap2 recipe](https://bioconda.github.io/recipes/minimap2/README.html)。

```bash
git clone https://github.com/xingguowu/IndelMarker.git
cd IndelMarker
conda env create -f envs/indelmarker.yaml
conda activate indelmarker
```

建议在创建前启用严格渠道优先级（Bioconda/conda-forge 的推荐配置）：`conda config --set channel_priority strict`。

然后在 `config.yaml` 中设置 `strategy: "1"`、`"2"` 或 `"3"`，填写相应输入路径，再运行：

```bash
snakemake --cores 16 --use-conda
```

策略 1：两个染色体级组装的 `minimap2 → SyRI → InDel 长度过滤 → Primer3`。结果为 `results/strategy1/filtered_indels.vcf` 及 `results/primer3/primer3_output.txt`。

策略 2：两个亲本 PE reads 的 `fastp → bwa → GATK HaplotypeCaller/联合分型 → 亲本基因型和长度过滤 → Primer3`。结果为 `results/strategy2/contrasting_indels.vcf` 及 Primer3 输出。样本固定为配置中的 `p1` 与 `p2`；模式 A 保留 0/0 对 1/1，模式 B 也保留杂合对 1/1。

策略 3：按 B73 注释在指定区间找基因，取最长 CDS，用 exonerate 在所选基因组中定位，提取带侧翼序列后与 CDS 用 MAFFT 比对。`annotation` 同时支持标准 GFF3 和提供的 `B73v5.gene.bed`。它不调用 Primer3；结果包括 `results/strategy3/alignment/pos.tsv`、`manifest.tsv` 与每个基因的 `.fa` 比对文件。

首次运行策略 1/2 前请确保参考 fasta 已可被 `samtools faidx`；策略 2 会自动生成 faidx、BWA 和 GATK 字典索引。工具名或绝对路径可在 `tools` 段覆盖。

安装完成后，先做环境和工作流解析检查：

```bash
for program in snakemake samtools minimap2 syri bcftools bwa gatk fastp primer3_core exonerate bedtools mafft; do command -v "$program"; done
snakemake --snakefile Snakefile --configfile config.yaml --dry-run --cores 1
```
