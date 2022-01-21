

            // Operator object //


            var operator_graphic = ['[{"accession": "EHN63435", "input_seq": "aggTCTTATGACATACGAaag", "lowest_perc_ident": "59.69", "num_seqs": 46, "consensus_score": 72.133, "motif": [{"base": "T", "score": 0.82}, {"base": "C", "score": 1.0}, {"base": "T", "score": 0.49}, {"base": "T", "score": 0.87}, {"base": "A", "score": 0.77}, {"base": "T", "score": 0.95}, {"base": "G", "score": 0.92}, {"base": "T", "score": 0.51}, {"base": "C", "score": 1.0}, {"base": "A", "score": 0.97}, {"base": "A", "score": 0.62}, {"base": "A", "score": 0.79}, {"base": "A", "score": 0.54}, {"base": "G", "score": 0.92}, {"base": "A", "score": 1.0}]}]']


            //Deals with inputting list of consensus motif data with different percent identity cutoffs

            for (i in range(0, operator_graphic.length)){
                var operator_data = JSON.parse(operator_graphic[i]);
                for (j in range(0, operator_data.length)){
                    renderOperator(operator_data[j]);
                }
            }

        