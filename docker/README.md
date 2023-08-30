# EnTAP Docker image Instructions
The following provides instructions for building and using the EnTAP docker image. 

## Build and Publish the Docker Image
Within this directory, run the `make` command for the version of EnTAP that you want to be build.
```bash
make v0.10.7-beta
```

If this is the official build of the EnTAP image to be shared on DockerHub then it can published to with the following command (only do this for official releases).
```bash
make push v0.10.7-beta
```

## Build a Docker Image with Local Development Code
To Dockerize the the current local development copy of your code you must run the make command in the root EnTAP directory with the following command:
```bash
make -f ./docker/Makefile local
```

## Setup output directory
EnTAP requires a one-time configuration step that processes all of the input files and downloads the EnTAP specific files. These resulting files will go into an `entap_outfiles` directory that must live outside of the Docker image because the files are too big for inclusion in the image.  

Create the output file directory so that it is not owned by root when Docker runs:
```bash
mkdir entap_outfiles
```

## Docker
### Test EnTAP
The following command will run the test data that comes with EnTAP.  Therefore, the `test_data` folder must also be mounted into the image.  The Docker image of EnTAP comes with a built-in config file for testing in the `/entap_config` directory of the image.

First, configure EnTAP with the test data.

```bash
docker run \
  --user "$(id -u)" \
  -v ${PWD}/test_data:/test_data \
  -v ${PWD}/entap_outfiles:/entap_outfiles \
  entap/entap:v0.10.7-beta \
  EnTAP --config \
    -t 8 \
    -d /test_data/swiss_prot_test.fasta \
    --out-dir /entap_outfiles \
    --ini /entap_config/entap_config.ini
```

Next run EnTAP
```bash
docker run \
  --user "$(id -u)" \
  -v ${PWD}/test_data:/test_data \
  -v ${PWD}/entap_outfiles:/entap_outfiles \
  entap/entap:v0.10.7-beta \
  EnTAP --runP \
    -i /test_data/trinity.fnn \
    -d /entap_outfiles/bin/swiss_prot_test.dmnd \
    --ini /entap_config/entap_config.ini
```

### Run EnTap with obtained data

#### Obtain Data Files
For these instructions we will download files into a folder named `input_files` which can be housed anywhere. For this tutorial we will place them in `/local/data/EnTAP/input_files`.  This where all EnTAP supported FASTA files should go.
```bash
mkdir -p /local/data/EnTAP/input_files
cd /local/data/EnTAP/input_files
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
```
Now that the image is constructed and we have the files downloaded, we will use the EnTAP docker image from the same directory where the `input_files` directory is housed.

Next we need a FASTA file containing transcripts. We can use the simple Tripal sweet orange example data:
```bash
mkdir entap_input
cd entap_input
wget http://tripal.info/sites/default/files/Citrus_sinensis-orange1.1g015632m.g.fasta
```

This data set requires the `-no_refine_starts` flag is set for TransDecoder. So we will provide a different INI file for EnTAP.  
```bash
cp ./docker/EnTAP-v0.10.7-beta/entap_config.ini ./entap_input/
```
Next edit the `entap_input/entap_config.ini` file and set the `transdecoder-no-refine-starts` setting to `true`

#### Run EnTAP
First, run the configuration step using `docker run`. The `input_files` and `entap_outfiles` directories are mounted within the image using the `-v` flag.  The image is setup such that EnTAP is configured to use these directories mounted directly to the root folder.

```bash
docker run \
  --user "$(id -u)" \
  -v /local/data/EnTAP/input_files:/input_files \
  -v ${PWD}/entap_input:/entap_input \
  -v ${PWD}/entap_outfiles:/entap_outfiles \
  entap/entap:v0.10.7-beta \
  EnTAP --config -t 8 \
    -d /input_files/uniprot_sprot.fasta \
    --out-dir /entap_outfiles \
    --ini /entap_input/entap_config.ini
```

Next run EnTAP on the data
```bash
docker run \
  --user "$(id -u)" \
  -v ${PWD}/entap_input:/entap_input \
  -v ${PWD}/entap_outfiles:/entap_outfiles \
  entap/entap:v0.10.7-beta \
  EnTAP --runP \
    -i /entap_input/Citrus_sinensis-orange1.1g015632m.g.fasta \
    -d /entap_outfiles/bin/uniprot_sprot.dmnd \
    --ini /entap_input/entap_config.ini
```

### Using EnTAP in Interactive Mode
You can work directly with EnTAP in interactive mode by running the docker image in "shell" mode. Once running, you can execute EnTAP commands.  You will be logged in as the `entap` user.

```bash
docker run \
  --user "$(id -u)" \
  -v /local/data/EnTAP/input_files:/input_files \
  -v ${PWD}/test_data:/test_data \
  -v ${PWD}/entap_outfiles:/entap_outfiles \
  -it \
  entap/entap:v0.10.7-beta \
  /bin/bash
```

## Singularity

```bash
singularity pull docker://entap/entap:v0.10.7-beta
```

### configuration

```bash
singularity exec \
  -B ${PWD}/entap_outfiles:/entap_outfiles,/local/data/EnTAP/input_files:/input_files \
  docker://entap/entap:v0.10.7-beta \
  EnTAP --config \
    -t 8 \
    -d /input_files/uniprot_sprot.fasta \
    --out-dir /entap_outfiles \
    --ini /entap_config/entap_config.ini

```

### Using EnTAP in Interactive Mode
```bash
singularity shell \
  -B /local/data/EnTAP/input_files:/input_files,${PWD}/entap_outfiles:/entap_outfiles \
  docker://entap/entap:v0.10.7-beta
```
