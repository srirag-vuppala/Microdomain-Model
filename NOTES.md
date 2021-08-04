Notes for the microdomain model
---

1. For the first pass, we'll consider epsilon = 0 so we can remove H from our Ionic term
2. How do we compute gap junctional resistance?
3. I think generate_E_matrix correctly generates E compared to create_E
    - Comparing to the notes we've made we defined C = [E 0] but the generate_E_matrix essentially makes the C matrix.
    - Yeah my bad create_E and create_C look good, I'm down to ditch the generate matrices functions since they
   do the same stuff as the create matrices functions
4. I think we still need to incorporate the gap junctional term as a diffusion term?
5. Why do we add constant to first and last entries of L3 matrix?
6. gin = 666 mS/cm^2, to non-dimensionalize this it looks like we need to muliplty by Rm, which is defined on p. 10 of weston's thesis ** ask dr. lin **
7. TO-DO: 
      - account for boundary cases in L2 and L3
      - change adding constant to multiplying at end points

