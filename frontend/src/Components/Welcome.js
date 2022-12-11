
import { Grid, Typography } from '@mui/material';

export default function Welcome() {

    return (


        <Grid
        container
        spacing={0}
        direction="column"
        alignItems="center" >

            <Typography
                sx={{ fontSize: { xs: 60, md: 90 }, mb: 1, mt: '5%' }}
                component="div" >
                GroovIO
            </Typography>

        </Grid>

    );
}