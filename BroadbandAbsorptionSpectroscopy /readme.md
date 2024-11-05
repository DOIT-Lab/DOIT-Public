Supporting code for:
Giles Blaney, Phillip Curtsmith, Angelo Sassaroli, Cristianne Fernandez, and Sergio Fantini, "Broadband absorption spectroscopy of heterogeneous biological tissue," Appl. Opt. 60, 7552-7562 (2021) 
https://doi.org/10.1364/AO.431013

# Data structure description for TissDataAPOPAI2021.mat

Data from 3 subjects stored in structs (sub1, sub2, and sub3),
organization of subX follows.

subX: Struct containing 9 fields

-   folderNames: Cell array specifying raw data location in DOIT Box
    > file structure {for internal use}

-   tissNames: Cell array specifying tissue names in data raw files {for
    > internal use}

-   plotTissNames: Cell array of tissue names

-   Abdomen / Breast / Forearm / Arm / Calf / Thigh: Struct arrays (over
    > repeated measurements) containing data from the tissue specified
    > by the struct name. These structs have the following fields

    -   FD: Struct containing frequency-domain data with the following
        > fields

        -   lambda: \[Number of wavelengths x 1\] array of wavelength.
            > (nm)

        -   rho: \[1 x Number of distances\] array of distances in the
            > dual-slope probe in the order \[1A, 1B, 2B, 2A\]. (mm)

        -   II: \[Number of wavelengths x Number of distances\] array of
            > frequency-domain amplitude (i.e., intensity) data.
            > (arb./mm\^2)

        -   PP: \[Number of wavelengths x Number of distances\] array of
            > frequency-domain phase data. (rad)

        -   II\_err: \[Number of wavelengths x Number of distances\]
            > array of error in frequency-domain amplitude (i.e.,
            > intensity) data. (arb./mm\^2)

        -   PP\_err: \[Number of wavelengths x Number of distances\]
            > array of error in frequency-domain phase data. (rad)

    -   CW: Struct containing continuous-wave data with the following
        > fields

        -   lambda: \[Number of wavelengths x 1\] array of wavelength.
            > (nm)

        -   rho: \[1 x Number of distances\] array of distances in the
            > dual-slope probe in the order \[1A, 1B, 2B, 2A\]. (mm)

        -   II: \[Number of wavelengths x Number of distances\] array of
            > intensity data. (arb./mm\^2)

        -   II\_err: \[Number of wavelengths x Number of distances\]
            > array of error in intensity data. (arb./mm\^2)