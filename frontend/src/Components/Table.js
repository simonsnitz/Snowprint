import { useState, useEffect } from 'react';
import { Box } from "@mui/material";
import { DataGrid } from '@mui/x-data-grid';
import { Link, Routes, Route } from 'react-router-dom';

import PleaseSelect from './PleaseSelect.js';
import Sensor from './Sensor.js';


export default function Table() {

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

        <Routes>
            <Route path="/" element={<PleaseSelect/>} />
            {sensorRouteList}
        </Routes>

        </Box>

    );
}