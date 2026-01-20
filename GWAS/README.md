GWAS流程

输入：vcf文件、phenotype文件

1. 1_GWAS.sh 文件中修改输入，运行
2. 将phenotype导入到生成的.fam文件中，推荐excel vlookup函数
3. 2_GWAS.sh 文件中修改输入，运行
4. 查看生成的 gwas_PCA.eigenval 文件，以确定该使用的PC数（awk提取重定向到新文件即可）
5. 运行 gemma
6. 可视化通过R脚本进行
