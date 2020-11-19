# EnTAP Docker image Instructions
The following provides instructions for building and using the EnTAP docker image.

## Build and Publish the Docker Image
Within this directory, use the docker command to build the image.  This command will create an EnTAP docker image that has all of the necessary prerequiesites.
```bash
docker build -t entap/entap:v0.10.7 .
```

If this is the official build of the EnTAP image then it can published to DockerHub with the following command (only do this for official releases).
```bash
docker push entap/entap:v0.10.7
```

## Obtain Data Files
For these instructions we will download files into a folder named `input_files` which can be housed anywhere. For this tutorial we will place them in `/local/data/EnTAP/input_files`.  This where all EnTAP supported FASTA files should go.
```bash
mkdir -p /local/data/EnTAP/input_files
cd /local/data/EnTAP/input_files
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
```
Now that the image is constructed and we have the files downloaded, we will use the EnTAP docker image from the same directory where the `input_files` directory is housed.

## Setup output directory
EnTAP requires a one-time configuration step that processes all of the input files and downloads the EnTAP specific files. These resulting files will go into an `entap_outfiles` directory that must live outside of the Docker image because the files are too big for inclusion in the image.  

Create the output file directory so that it is not owned by root when Docker runs:
```bash
mkdir entap_outfiles
```

## Docker
### Configure EnTAP
Next, run the configuration step using `docker run`. The `input_files` and `entap_outfiles` directories are mounted within the image using the `-v` flag.  The image is setup such that EnTAP is configured to use these directories mounted directly to the root folder.

```bash
docker run \
  --user "$(id -u)" \
  -v /local/data/EnTAP/input_files:/input_files \
  -v ${PWD}/test_data:/test_data \
  -v ${PWD}/entap_outfiles:/entap_outfiles \
  entap/entap:v0.10.7 \
  EnTAP --config -d /input_files/uniprot_sprot.fasta --ini /home/entap/entap_config.ini
```

### Test EnTAP
The following command will run the test data that comes with EnTAP.  Therefore, the `test_data` folder must also be mounted into the image.  

```bash
docker run \
  --user "$(id -u)" \
  -v ${PWD}/test_data:/test_data \
  -v ${PWD}/entap_outfiles:/entap_outfiles \
  entap/entap:v0.10.7 \
  EnTAP --config -d /test_data/swiss_prot_test.fasta --out-dir /test_data --ini /home/entap/entap_config.ini
```

```bash
docker run \
  --user "$(id -u)" \
  -v ${PWD}/test_data:/test_data \
  -v ${PWD}/entap_outfiles:/entap_outfiles \
  entap/entap:v0.10.7 \
  EnTAP --runP -i /test_data/trinity.fnn -d /test_data/bin/swiss_prot_test.dmnd --ini /home/entap/entap_config.ini
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
  entap/entap:v0.10.7 \
  /bin/bash
```

## Singularity

```bash
singularity pull docker://entap/entap:v0.10.7
```

## Using EnTAP in Interactive Mode
```bash
singularity shell \
  -B /local/data/EnTAP/input_files:/input_files,${PWD}/entap_outfiles:/entap_outfiles \
  docker://entap/entap:v0.10.7
```
