import {Box} from "@mui/material";
import Header from "./Components/Header.js";
import AllSensorsTable from "./Components/AllSensorsTable.js";

import { Routes, Route } from 'react-router-dom';

function App() {
  return (
    <Box 
      sx={{backgroundColor:"white"}}>

        <Header/>

          {/* Probably unecessary to use Routes here */}
        <Routes>
          <Route path="/*" element={<AllSensorsTable />} />
        </Routes>

    </Box>
  );
}

export default App;
