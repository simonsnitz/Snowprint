import {Box} from "@mui/material";
import Welcome from "./Components/Welcome.js";
import Table from "./Components/Table.js";


function App() {
  return (
    <Box 
      sx={{backgroundColor:"#edf4ff"}}>

        <Welcome/>
        <Table
        />

    </Box>
  );
}

export default App;
