import { useState, useEffect } from 'react';
import { Box } from "@mui/material";
import { DataGrid } from '@mui/x-data-grid';
import Logo from './Logo';


export default function Table() {

const [data, setData] = useState([]);
const [rows, setRows] = useState([]);

//const [protein, setProtein] = useState([]);

useEffect(()=>{

    fetch('data.json'
    , {
        headers : { 
            'Content-Type': 'application/json',
            'Accept': 'application/json',
        }
    }
    )
    .then((res) => res.json())
    .then((data) => {
        console.log(data);
        setData(data);
    });
    
}, [])

useEffect(()=>{

    const rowsToAdd = [];

    if (typeof data !== 'undefined') {
        let counter = 0;
        for (var entry in data) {
            var entry = {
                id: counter,
                accession: data[entry].accession,
                score: data[entry].score,
                sequencesAligned: data[entry].sequencesAligned,
                organism: data[entry].organism
            };
            rowsToAdd.push(entry);
        
            // Add code here to add Routes

        counter += 1;
    }


        setRows(rowsToAdd);
    }
  },[data])



const columns = [
    {   
        field: 'id', 
        headerName: 'ID', 
        width: 70 },
    { 
        field: 'accession', 
        headerName: 'Accession', 
        width: 200 },
    {
        field: 'score',
        headerName: 'Score',
        type: 'number',
        width: 70,
    },
    {
        field: 'sequencesAligned',
        headerName: 'Sequences Aligned',
        width: 200,
    },
    { 
        field: 'organism', 
        headerName: 'Organism', 
        width: 200 },
  ];

// const rows = [
//     { id: 1, accession: data.accession, score: data.score, sequencesAligned: data.sequencesAligned, organism: data.organism }
// ];




    return (


        <Box 
            sx={{ height: 400, width: '90%', backgroundColor: "white", margin: "auto" }}
        >

            <DataGrid
                    rows={rows}
                    columns={columns}
                    rowsPerPageOptions={[10]}
                    density="compact"
                    autoPageSize
            />

            <Logo
                data="PROTEIN"
            />

        </Box>

    );
}