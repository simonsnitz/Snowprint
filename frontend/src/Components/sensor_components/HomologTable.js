import { useState, useEffect } from 'react';
import { Box } from "@mui/material";
import { DataGrid } from '@mui/x-data-grid';


export default function HomologTable(props) {


const [rows, setRows] = useState([]);
const data = props.data;

useEffect(() => {

    const rowsToAdd = [];

    if (typeof props.data !== 'undefined') {
        let counter = 0;
        for (var entry in data) {
            var protein = {
                id: counter,
                accession: data[entry].accession,
                identity: data[entry].identity,
                coverage: data[entry].coverage,
                predicted_operator: data[entry].predicted_operator,
                align_score: data[entry].align_score,
                organism: data[entry].organism
            };
            rowsToAdd.push(protein);

            counter += 1;
        }
    
    
            setRows(rowsToAdd);
        }
      },[props.data, data])



const columns = [
    {   
        field: 'id', 
        headerName: 'ID', 
        width: 70 },
    { 
        field: 'accession', 
        headerName: 'Accession', 
        width: 140,
        renderCell: (params) => (
            <a href={"https://www.ncbi.nlm.nih.gov/protein/"+params.value} 
                target="__blank">
                {params.value}
            </a>
        ) 
    },
{ 
        field: 'identity', 
        headerName: 'Protein identity', 
        type: 'number',
        width: 110 },
{ 
        field: 'coverage', 
        headerName: 'Protein coverage', 
        type: 'number',
        width: 120 },
    { 
        field: 'predicted_operator', 
        headerName: 'Predicted operator', 
        width: 400 },
    {
        field: 'align_score',
        headerName: 'Alignment score',
        type: 'number',
        width: 120,
    },
    { 
        field: 'organism', 
        headerName: 'Organism', 
        width: 200 },
  ];




    return (


        <Box
            sx={{ height: 600, width: '90%', backgroundColor: "white", margin: "auto" }}
        >

            <DataGrid
            sx={{ mb:50 }}
                    rows={rows}
                    columns={columns}
                    //rowsPerPageOptions={[100]}
                    density="compact"
                    //autoPageSize
            />

        </Box>

    );
}