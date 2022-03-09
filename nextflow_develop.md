nextflow流程开发


1.	nextflow流程目标用户：熟悉shell常用命令的研究人员？

2.	开发流程需由工具开发相关人员提供：
(1)	每个工具步骤的输入输出，参数详解；

(2)	源码：github; conda; pypi for python packages; r-cran/Bioconductor for R packages

3.	环境配置
(1)	conda: 如可用environment.yaml管理依赖

conda env create -f environment.yml
install guppy

…
(2)	docker: 利用docker配置依赖和自安装包
…

4.	nextflow流程开发
(1)	环境依赖及其它配置，与第3步结合

(2)	流程开发，结合第二步合理规划流程步骤，编写nextflow语句。最后可形成DAG图？


5.	测试


（开发过程中的困难/想法可一一列出）




需熟悉shell












