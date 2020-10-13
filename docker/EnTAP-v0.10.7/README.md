# EnTAP Docker image Instructions
The following provides instructions for building the EnTAP docker image.

## Step 1. Build the image
Within this directory, use the docker command to build the image.  This command will create an EnTAP docker image that has all of the necessary prerequiesites.
```bash
docker build -t entap-docker .
```

## Step 2. Obtain data files.
For these instructions we will download files into a folder named `input_files` which is housed one level above this directory.  This where all EnTAP supported FASTA files should go.
```bash
mkdir ../input_files
cd ../input_files
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
```
Now that the image is constructed and we have the files downloaded, we will use the EnTAP docker image from the same directory where the `input_files` directory is housed.

## Step 3. Configure EnTAP
EnTAP requires a one-time configuration step that processes all of the input files and downloads the EnTAP specific files. These resulting files will go into an `entap_outfiles` directory that must live outside of the Docker image because the files are too big for inclusion in the image.  

First, create the output file directory so that it is not owned by root when Docker runs:
```bash
mkdir entap_outfiles
```

Second, run the configuration step using `docker run`. The `input_files` and `entap_outfiles` directories are mounted within the image using the `-v` flag.  The image is setup such that EnTAP is configured to use these directories mounted directly to the root folder.

```bash
docker run --user "$(id -u):$(id -g)" -v ${PWD}/input_files:/input_files -v ${PWD}/entap_outfiles:/entap_outfiles entap-docker EnTAP --config -d /input_files/uniprot_sprot.fasta -t 8
```

## Step 4. Test EnTAP
The following command will run the test data that comes with EnTAP.  Therefore, the `test_data` folder must also be mounted into the image.  
```bash
docker run --user "$(id -u):$(id -g)" -v ${PWD}/test_data:/test_data -v ${PWD}/entap_outfiles:/entap_outfiles entap-docker EnTAP --runP -i /test_data/trinity.fnn
```

## Step 5. (Optional) Use EnTAP in Interactive Mode
You can work directly with EnTAP in interactive mode by running the docker image in "shell" mode. Once running, you can execute EnTAP commands.  You will be logged in as the `entap` user.

```bash
docker run --user "$(id -u)" -v ${PWD}/input_files:/input_files -v ${PWD}/test_data:/test_data -v ${PWD}/entap_outfiles:/entap_outfiles -it entap-docker /bin/bash
```
