import {Box} from "@mui/material";
import Welcome from "./Components/Welcome.js";
import Table from "./Components/Table.js";

import { Routes, Route } from 'react-router-dom';

function App() {
  return (
    <Box 
      sx={{backgroundColor:"#edf4ff"}}>

        <Welcome/>

          {/* Probably unecessary to use Routes here */}
        <Routes>
          <Route path="/*" element={<Table />} />
        </Routes>

    </Box>
  );
}

export default App;
