## nextflow流程开发

以DeepSignal流程为例

<p>&nbsp;</p>

### 1. nextflow流程目标用户：熟悉shell常用命令

<p>&nbsp;</p>

### 2. 开发流程需由工具开发相关人员提供：

#### (1) 每个工具步骤的输入输出，参数详解

DeepSignal步骤：

# ![deepsignal-pipeline](images/pipeline_dag_2022-03-09_17-04-28.png)

# ![deepsignal-quickstart](images/deepsignal-quickstart.png)
参考：[https://github.com/bioinfomaticsCSU/deepsignal#quick-start](https://github.com/bioinfomaticsCSU/deepsignal#quick-start)

#### (2) 源码
  * __github__

  * __conda__/bioconda（推荐，方便环境配置） 

  * __pypi__ (python packages)

# ![deepsignal-pypi](images/deepsignal-pypi.png)

  * __r-cran/Bioconductor__ (R packages) 

<p>&nbsp;</p>

### 3. nextflow流程开发

#### (1) 环境依赖配置

  - conda

用environment.yml管理依赖。

示例：
# ![longmethyl-environment](images/longmethyl-environment.png)

参考：[longmethyl/environment.yml](../environment.yml)

  - docker（推荐）

利用docker配置依赖和自安装包，编写Dockerfile。

```sh
cd longmethyl/
# run docker build
docker build .
# or add tags
docker build -t nipengcsu/longmethyl:0.1 .
# push to dockerhub [create repo in docker hub first, and login required]
docker push nipengcsu/longmethyl:0.1
```

示例：
# ![longmethyl-dockerfile](images/longmethyl-dockerfile.png)

参考：[longmethyl/Dockerfile](../Dockerfile)

#### (2) 流程开发，编写nextflow语句

从[main.nf](../main.nf)开始编写nextflow流程，配置参数在[nextflow.config](../nextflow.config).

  - process举例：

# ![nextflow-process](images/nextflow-process.png)

  - workflow举例：

# ![nextflow-workflow](images/nextflow-workflow.png)

<p>&nbsp;</p>

### 4. 测试

#### (1) 安装nextflow
#### (2) 下载longmethyl流程源码
#### (3) 准备数据，运行程序

deepsignal流程需准备：fast5 files, genome fasta, model.ckpt

# ![longmethyl-test](images/longmethyl-test.png)

详见[README.md](../README.md)。
