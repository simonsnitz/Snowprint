import { Box, Grid, Typography } from '@mui/material';

import Logo from './Logo.js';

export default function Sensor(props) {

    return (

        <Box>

            <Grid container spacing={3} sx={{ mt: 2 }} justifyContent="center">

                {/* Alias  */}
                <Grid item xs={12}>
                    <Typography
                    component="div"
                    gutterBottom
                    sx={{ fontSize: { xs: 30, sm: 50 }, textAlign: 'center' }}
                    >
                    {props.accession}
                    </Typography>
                </Grid>

                {/* Logo */}
                <Grid item xs={12}>
                    <Logo 
                        data = {props.data}
                        />
                </Grid>

            </Grid>

        </Box>


    );
}