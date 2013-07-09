     #include <iostream>
     #include <stdlib.h>
     #include <stdio.h>
     #include <netcdf.h>
     
     /* Handle errors by printing an error message and exiting with a
      * non-zero status. */
     #define ERRCODE 2
     #define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
     
     int
     main(int argc, char *argv[])
     {
        using namespace std;
        /* This will be the netCDF ID for the file and data variable. */
        int ncid, varid, dimid;
        int ndims_in, nvars_in, ngatts_in, unlimdimid_in;

        /* Loop indexes, and error handling. */
        int x, y, retval ;
        size_t x_val, y_val ;
     
        if (argc <= 3)
        {
           std::cout << "Usage: " << argv[0] << " <Filename> <x_dim> <y_dim>" << std::endl;
           exit(1);
        }

        char *pFilename = argv[1];
        int  x_true = atoi(argv[2]); 
        int  y_true = atoi(argv[3]); 
 
	cout << "Processing: " << pFilename << " x " << x_true << " y " << y_true << endl;

        /* Open the file. NC_NOWRITE tells netCDF we want read-only access
         * to the file.*/
        if ((retval = nc_open(pFilename, NC_NOWRITE, &ncid)))
           ERR(retval);
     
        /* There are a number of inquiry functions in netCDF which can be
           used to learn about an unknown netCDF file. NC_INQ tells how
           many netCDF variables, dimensions, and global attributes are in
           the file; also the dimension id of the unlimited dimension, if
           there is one. */
        if ((retval = nc_inq(ncid, &ndims_in, &nvars_in, &ngatts_in, &unlimdimid_in)))
           ERR(retval);

        // cout << "ndims " << ndims_in << "nvars " << nvars_in << endl;

        /* Get the varid of the data variable, based on its name. */

        if ((retval = nc_inq_dimid(ncid, "x", &dimid)))
           ERR(retval);

        if ((retval = nc_inq_dimlen(ncid, dimid, &x_val )))
           ERR(retval);

        if ((retval = nc_inq_dimid(ncid, "y", &dimid)))
           ERR(retval);

        if ((retval = nc_inq_dimlen(ncid, dimid, &y_val )))
           ERR(retval);

        // cout << "x " << x_val << "y " << y_val << endl;

        /* Close the file, freeing all resources. */
        if ((retval = nc_close(ncid)))
           ERR(retval);
     
        retval =  !((x_true == x_val)&&(y_true==y_val));
        cout << "retval " << retval << endl;

        return retval;
     }

