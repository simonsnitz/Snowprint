import { Box, Grid } from '@mui/material';
import { DNALogo } from 'logojs-react';


export default function Snowprint(props) {


// Collect all operators. Assert same length and only ATCG characters
  const operators = [];
  const op_length = props.data[0]["predicted_operator"].length;

  for (const entry of props.data) {
    var op = entry["predicted_operator"].toUpperCase();
    if (op.length === op_length && op.match(/^[ATCG]*$/)){
      operators.push(op);
      }
    }



// Create the base pair frequency matrix
  const num_ops = operators.length
  const matrix = []

    // define the position
  for (let i = 0; i < operators[0].length; i++){

    var base = [0,0,0,0];
      // loop through each operator
      for (const op of operators){
        if (op[i] === "A"){
          base[0] += 1/num_ops
        }
        else if (op[i] === "C"){
          base[1] += 1/num_ops
        }
        else if (op[i] === "G"){
          base[2] += 1/num_ops
        }
        else if (op[i] === "T"){
          base[3] += 1/num_ops
        }
      }
      base[0] = Math.round((base[0]+ Number.EPSILON) * 100) / 100;
      base[1] = Math.round((base[1]+ Number.EPSILON) * 100) / 100;
      base[2] = Math.round((base[2]+ Number.EPSILON) * 100) / 100;
      base[3] = Math.round((base[3]+ Number.EPSILON) * 100) / 100;

      matrix.push(base)
  }

console.log(operators);



  return (
    <Box>
      <Grid container spacing={4} columns={12} justifyContent="center">
        <Grid item xs={12}>

            <DNALogo ppm={matrix} overflow="scroll"
                height="200px"/>

        </Grid>
      </Grid>
    </Box>
  );
}