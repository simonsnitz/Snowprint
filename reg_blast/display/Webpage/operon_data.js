

            // Operon object //

            var operon_graphic = [{'meta_data': '{"accession": "WP_003080639", "identity": 100.0, "organism": "placeholder"}', 'operon_data': '{"0": {"length": "10.0%", "description": "acyl-CoA dehydrogenase family protein", "color": "#ff4242", "direction": "-", "spacer": "0.0%", "link": "WP_003080627.1"}, "1": {"length": "17.0%", "description": "biotin/lipoyl-binding protein", "color": "black", "direction": "-", "spacer": "0.0%", "link": "WP_003080628.1"}, "2": {"length": "7.0%", "description": "enoyl-CoA hydratase-related protein", "color": "#ff4242", "direction": "-", "spacer": "0.0%", "link": "WP_003080630.1"}, "3": {"length": "10.0%", "description": "acyl-CoA dehydrogenase family protein", "color": "#ff4242", "direction": "-", "spacer": "0.0%", "link": "WP_003080632.1"}, "4": {"length": "15.0%", "description": "acetyl-CoA carboxylase carboxyltransferase subunit", "color": "#ff4242", "direction": "-", "spacer": "0.0%", "link": "WP_003080634.1"}, "5": {"length": "8.0%", "description": "SDR family oxidoreductase", "color": "#ff4242", "direction": "-", "spacer": "0.0%", "link": "WP_003080637.1"}, "6": {"length": "6.0%", "description": "TetR/AcrR family transcriptional regulator", "color": "blue", "direction": "-", "spacer": "1.0%", "link": "WP_003080639.1"}, "7": {"length": "9.0%", "description": "AraC family transcriptional regulator", "color": "#6eff42", "direction": "-", "spacer": "1.0%", "link": "WP_003080641.1"}, "8": {"length": "10.0%", "description": "lipid-transfer protein", "color": "black", "direction": "+", "spacer": "0.0%", "link": "WP_003080643.1"}}'}]

            var len = operon_graphic.length
            for (i in range(0, len)){
                var metaData = JSON.parse(operon_graphic[i]["meta_data"]);
                var operon = JSON.parse(operon_graphic[i]["operon_data"]);
                renderOperon(metaData, operon);
            }
                    