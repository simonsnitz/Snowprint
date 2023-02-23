[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/simonsnitz/Snowprint/)
[![Issues](https://img.shields.io/github/issues/althonos/pyfamsa.svg?style=flat-square&maxAge=600)](https://github.com/simonsnitz/Snowprint/issues)


# Snowprint [![Stars](https://img.shields.io/github/stars/simonsnitz/Snowprint.svg?style=social&maxAge=3600&label=Star)](https://github.com/simonsnitz/Snowprint/stargazers)


## 🗺️ Overview

Snowprint is a bioinformatic tool used to predict the DNA-binding sequence of transcription factors. As opposed to DNAse footprinting, it requires no experimentation.



## 🔧 Installing

1. Clone the repo
```console
$ git clone https://github.com/simonsnitz/Snowprint.git
```

2. Setup the virtual environment
```console
$ python -m venv env
$ source env/bin/activate
$ pip install -r requirements.txt
```

3. Install blast

*On Mac*
```console
$ brew install blast
```
*On Linux*
```console
$ sudo apt-get install ncbi-blast+
```

4. Setup the frontend interface
#### Download node
https://nodejs.org/en/download/

#### Install node packages
```console
$ cd Snowprint/frontend
$ npm install
```
#### Launch the frontend interface
```console
$ npm start
```
*This will open the frontend interface in your browser*

## Create a Snowprint prediction
```console
$ cd Snowprint
$ python main.py {your RefSeq accession goes here}
```



## Caching
- A SQLite database is used to cache intermediate data collected while creating a prediction to increase prediction speed.
- The database is located in `Snowprint/Cache/Snowprint.db`
- The database schema can be found in `Snowprint/Cache/models.py`


## Display data
- All data for frontend display is located in `Snowprint/frontend/public/data.json`