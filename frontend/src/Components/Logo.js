import { useState, useEffect } from 'react';
import { Grid, Typography } from '@mui/material';

export default function Logo(props) {

    return (

        <Grid
        container
        spacing={0}
        direction="column"
        alignItems="center" >

            <Typography
                sx={{ fontSize: { xs: 20}, mb: 1, mt: '5%' }}
                component="div" >
                
                {props.data}
                
            </Typography>

        </Grid>

    );
}