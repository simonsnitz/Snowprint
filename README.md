# GroovIO
Predicts the operator sequence for prokaryotic transcrition regulators.

## Usage
#### Perform analysis on a regulator

1. Clone the repo
```
$ git clone https://github.com/simonsnitz/GroovIO.git
```

2. Activate the virtual environment
```
$ source env/bin/activate
```

3. Install dependencies
```
$ pip -r requirements.txt
```

4. Input your regulator accession ID(s) in the GroovIO/cache/reg_ids/reg_ids.txt file

5. Run the prediction algorithm on your regulators. These will be added to the SQL database
```
$ python main.py
```

## Dependencies
- python3-venv
- python3-pip
- biopython
- requests
- sqlalchemy


## Notes
This is a tool used by [The Biosensor Database](https://gbiosensors.com).

## License

Work in the `src` directory is licensed under the MIT License, found in the LICENSE.txt file

