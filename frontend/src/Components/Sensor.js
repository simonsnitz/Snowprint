import { Box, Grid, Typography, Paper } from '@mui/material';

import Snowprint from './sensor_components/Snowprint.js';
import HomologTable from './sensor_components/HomologTable.js';

export default function Sensor(props) {

    return (

        <Box>

            <Grid container spacing={3} sx={{ mt: 2 }} justifyContent="center">

                {/* Alias  */}
                <Grid item xs={12}>

                </Grid>

                {/* Logo & Metadata */}
                <Grid item xs={12} mb={5}>
                    <Snowprint data={props.data} />

                    <Paper elevation={5} sx={{ padding: 3, mb: 5 }}>
                        <Typography
                        component="div"
                        sx={{ fontSize: { xs: 12, sm: 16 }, overflowWrap: 'anywhere' }}
                        >
                        {props.intergenic}
                        </Typography>
                    </Paper>

                    <HomologTable data={props.data}/>
                </Grid>

            </Grid>

        </Box>


    );
}