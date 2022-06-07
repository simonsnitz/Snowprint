    # Create a list of representative protein accession IDs from a search term

    # Example
        # Inputs:
            # term:                                 TetR
            # amino acid length filter:             100-300
            # clustering percent identity:          50
            # min number of homologs per cluster:   20
        # Output:
            # text file with 9,577 accession IDs




    # 1. fetch fasta sequences using the search term
#query2db.py


    # 2. filter fasta sequences by amino acid length range
# need to use online galaxy tool for this at the moment


    # 3. cluster by X% identity
# command line argument using cd-hit


    # 4. filter clusters containing a minimum number of homologs
# filterClusters_minHomologs.py


