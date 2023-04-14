import { Box, Grid, Typography, Paper } from '@mui/material';

import Snowprint from './sensor_components/Snowprint.js';
import HomologTable from './sensor_components/HomologTable.js';

export default function Sensor(props) {

    return (

        <Box>

            <Grid container spacing={3} sx={{ mt: 5 }} justifyContent="center">

                {/* Alias  */}
                <Grid item xs={12}>
                <Typography
                        sx={{ fontSize: { xs: 18, md: 36 }, fontWeight: 300, textAlign: 'center' }}>
                        {props.accession}
                    </Typography>
                </Grid>

                {/* Logo & Metadata */}
                <Grid item xs={12}>

                    <Grid item mt={3} overflow="scroll">
                    <Typography
                        sx={{ fontSize: { xs: 18, md: 22 }, fontWeight: 300, textAlign: 'left', mb:2}}>
                        Snowprint
                    </Typography>
                    <Box justifyContent="center"  display="flex">
                    <Snowprint data={props.data} />
                    </Box>
                    </Grid>


                    <Grid item mt={10}>
                    <Typography
                        sx={{ fontSize: { xs: 18, md: 22 }, fontWeight: 300, textAlign: 'left', mb:1 }}>
                        Inter-operon region
                    </Typography>

                    <Paper elevation={5} sx={{ padding: 3, mb: 5 }}>
                        <Typography
                        component="div"
                        sx={{ fontSize: { xs: 12, sm: 16 }, overflowWrap: 'anywhere' }}
                        >
                        {props.intergenic}
                        </Typography>
                    </Paper>
                    </Grid>

                    <Grid item mt={10} mb={10}>
                    <Typography
                        sx={{ fontSize: { xs: 18, md: 22 }, fontWeight: 300, textAlign: 'left', mb:1 }}>
                        All homologs
                    </Typography>
                    <HomologTable data={props.data}/>
                    </Grid>
                </Grid>

            </Grid>

        </Box>


    );
}