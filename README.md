iSynthesis
=========
Forward synthesis planning system based on Monte-Carlo Tree Search


Using iSynthesis via docker
=========
Download docker [image](http://seafile.cimm.site/f/dbc828a4753e46d0987c/)

### Setup and run:

1. Loading docker image:

```
cat isynthesis.tar.gz | docker load
```

2. SDF file with target molecule (ex: bevantolol.sdf) should exist in your workdir. Run synthesis via python within docker container:

```
sudo docker run -it –rm -e reagents=1000 -e steps=3000 -e target=bevantolol -e cpu=4 -v $PWD:/mnt isynthesis:latest
```
Customize settings:

- reagents - number of building blocks to start from
- steps - number of MCTS iterations
- target - name of SDF file with target molecule
- cpu - numper of CPU to use

Every 100 iterations the results will be updated to a text file results.txt in the current folder. The process of searching for synthesis paths can be monitored in the log file in current directory.

Setup manually
=========
### Configure project via Poetry 

  ```
  pip install poetry
  poetry config virtualenvs.in-project true
  ```
  Create virtual environment and install dependencies:
  ```
  poetry shell
  poetry install
  ```
### Data
Firstly, download [building blocks](http://seafile.cimm.site/f/da3dba6f24694018bd40/) and

Extract bzip2 archive:
```
tar xjvf ZINCbb.tar.bz2
```

Secondly, download [reaction rules](https://seafile.cimm.site/f/891c595f7b7c422b9ae5/?dl=1) and

Extract tar.gz:
```
tar xvf rules.tar.gz 
```

To create your custom rules, follow the [tutorial-rules](https://github.com/Pandylandy/iSynthesis/blob/master/iSynthesis/tutorial/tutorial-rules.ipynb)

Then prepare CGRdb cartridge by [instructions](https://github.com/Pandylandy/CGRdb#readme)

To download building blocks from ZINC database, follow the [tutorial-zinc](https://github.com/Pandylandy/iSynthesis/blob/master/iSynthesis/tutorial/tutorial-zinc.ipynb)

To create database with downloaded building blocks, follow the [tutorial-database](https://github.com/Pandylandy/iSynthesis/blob/master/iSynthesis/tutorial/tutorial-database.ipynb)
