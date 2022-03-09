## nextflow流程开发

以DeepSignal流程为例

<p>&nbsp;</p>

### 1. nextflow流程目标用户：需熟悉shell常用命令

<p>&nbsp;</p>

### 2. 开发流程需由工具开发相关人员提供：

<p>&nbsp;</p>

#### (1) 每个工具步骤的输入输出，参数详解；

<p>&nbsp;</p>
#### (2) 源码：github; conda; pypi for python packages; r-cran/Bioconductor for R packages

<p>&nbsp;</p>

### 3. nextflow流程开发

<p>&nbsp;</p>
#### (1) 环境依赖配置

<p>&nbsp;</p>
##### conda: 可用environment.yml管理依赖

```sh
conda env create -f environment.yml
# install guppy also
```

<p>&nbsp;</p>
##### docker: 利用docker配置依赖和自安装包


<p>&nbsp;</p>
#### (2) 流程开发，结合第二步合理规划流程步骤，编写nextflow语句

<p>&nbsp;</p>
### 4. 测试













