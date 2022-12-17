import { Box, Grid, Typography, Paper } from '@mui/material';

import Logo from './Logo.js';
import MetadataTable from './MetadataTable.js';

export default function Sensor(props) {

    return (

        <Box>

            <Grid container spacing={3} sx={{ mt: 2 }} justifyContent="center">

                {/* Alias  */}
                <Grid item xs={12}>

                </Grid>

                {/* Logo & Metadata */}
                <Grid item xs={12} mb={5}>
                    <Logo data={props.data} />

                    <Paper elevation={5} sx={{ padding: 3, mb: 5 }}>
                        <Typography
                        component="div"
                        sx={{ fontSize: { xs: 12, sm: 16 }, overflowWrap: 'anywhere' }}
                        >
                        {props.intergenic}
                        </Typography>
                    </Paper>

                    <MetadataTable data={props.data}/>
                </Grid>

            </Grid>

        </Box>


    );
}