import { useState, useEffect } from 'react';
import { Box, Grid, Typography } from "@mui/material";
import { DataGrid } from '@mui/x-data-grid';
import { Link, Routes, Route } from 'react-router-dom';


import Sensor from './Sensor.js';


export default function AllSensorsTable() {

const [data, setData] = useState([]);
const [rows, setRows] = useState([]);
const [sensorRouteList, setSensorRouteList] = useState(null);

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

useEffect(() => {

    const rowsToAdd = [];
    const sensorRouteList = [];

    if (typeof data !== 'undefined') {
        let counter = 0;
        for (var entry in data) {
            var protein = {
                id: counter,
                accession: data[entry].accession,
                score: data[entry].score,
                sequencesAligned: data[entry].sequencesAligned,
                organism: data[entry].organism
            };
            rowsToAdd.push(protein);
        
           
            sensorRouteList.push(
                <Route
                  key={counter}
                  path={data[entry].accession}
                  element={
                    <Sensor
                        accession = {data[entry].accession} 
                        score = {data[entry].score}
                        sequencesAligned = {data[entry].sequencesAligned} 
                        organism = {data[entry].organism}
                        intergenic = {data[entry].intergenic}
                        data = {data[entry].data}
                    />
                  }
                />
              );


        counter += 1;
    }


        setRows(rowsToAdd);
        setSensorRouteList(sensorRouteList);
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
        width: 200,
        renderCell: (params) => (
            <Link to={params.value} >
                {params.value}
            </Link>
        ) },
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


  const selectionPrompt = () => {
    return (
      <Box>
        <Grid container spacing={4} columns={12} mt={5} justifyContent="center">
          <Grid item xs={10} mb={6}>
            <Typography
              sx={{ fontSize: { xs: 22, md: 24 }, textAlign: 'center' }}
            >
              Please select a protein
            </Typography>
          </Grid>
        </Grid>
      </Box>
    );
  };



    return (


        <Box 
            sx={{ height: 400, width: '90%', backgroundColor: "white", margin: "auto" }}
        >
            <Typography
              sx={{ fontSize: { xs: 18, md: 22 }, fontWeight: 300, textAlign: 'left', mb:2 }}
            >
              Processed regulators
            </Typography>
            <DataGrid
                    rows={rows}
                    columns={columns}
                    rowsPerPageOptions={[10]}
                    density="compact"
                    autoPageSize
            />

        <Routes>
            <Route path="/" element={selectionPrompt()} />
            {sensorRouteList}
        </Routes>

        </Box>

    );
}