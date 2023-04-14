
import { Grid, Typography, Box } from '@mui/material';


export default function Header() {

    return (


        <Grid
        container
        spacing={0}
        columns={12}
        mt={8}
        sx={{ mb:{xs:8, sm:0}}}
        alignItems="center" >


            <Grid item xs={0} sm={2} md={4}></Grid>

            <Grid item xs={4} sm={2} md={1}>
            <img src="Logo.png" alt="logo" width="100px" height="100px"/>
            </Grid>

            <Grid item xs={6} sm={3}>
            <Typography
                sx={{ fontSize: { xs: 50, sm: 70, md: 90 }, fontWeight: 300, mb: 5, mt: -3 }}
                display="inline" >
                Snowprint
            </Typography>
            </Grid>

        </Grid>

    );
}