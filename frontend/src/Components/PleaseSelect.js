import React from 'react';

import { Box, Grid, Typography } from '@mui/material';


export default function PleaseSelect() {
  return (
    <Box>
      <Grid container spacing={4} columns={12} mt={8} justifyContent="center">
        <Grid item xs={10} mb={6}>
          <Typography
            sx={{ fontSize: { xs: 22, md: 24 }, textAlign: 'center' }}
          >
            Please select a protein
          </Typography>
        </Grid>
      </Grid>
    </Box>
  );
}