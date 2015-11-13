/********************************************************************** 
** Copyright, 1998 - 2014, Intergraph Corporation. All rights reserved.
***********************************************************************/

/** @file dexample1.c */

/**
 * @page Examples
 * @section Decompression Decompression
 * @subsection dexample1 Decompression Example 1
 * @link dexample1.c @endlink <br>
 * 
 * This example demonstrates the <i><b>blocking</b></i> interface into the NCSEcw
 * library.  The application opens a view, reads the view, then reads another view.
 * 
 * This example uses the BIL read call; you could instead use the RGB call if
 * you want to the library to always return a straight RGB image regardless of source.
 * 
 * This example demonstrates the @link NCSOpenFileView @endlink, 
 * @link NCSCloseFileView @endlink, @link NCSGetViewFileInfo @endlink, 
 * @link NCSSetFileView @endlink and @link NCSReadViewLineBIL @endlink functions.
 */

#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include "NCSECWClient.h"
#include "NCSErrors.h"
#include "NCSString.h"
#include "NCSMisc.h"

#define MAX_WINDOW	333			// max size of a window view for testing
#define MAX_REGION_READS 10		// max number of reads of a file view for a client
#define RAND() ((double) rand() / (double) RAND_MAX)

int setup_report(int argc, char *argv[],
			char **p_p_input_ecw_filename);
int test_openview( char *szInputFilename, BOOLEAN bRandomReads, BOOLEAN bReportTime );

int main(int argc, char **argv)
{

	char	*szInputFilename;
	int		nError = 0;

	/*
	 * 	Initialize the library if we are linking statically
	 */
	NCSInit();
	
  	/*
	**	Find file names, and open the output
	*/
	if( setup_report(argc, argv,
		&szInputFilename) )
	{
		return( 1 );
	}

	nError = test_openview(szInputFilename, FALSE, TRUE);
	if( nError ) {
		printf("Openview test returned an error\n");
		NCSShutdown();        
		return(nError);
	}

	NCSShutdown();

	return(0);
}



int setup_report(int argc, char *argv[],
			char **ppInputFilename)
{
	if (argc != 2) {
	  printf("Usage: %s <input filename.ecw>\n", argv[0]);
	  return(1);
	}

	*ppInputFilename = argv[1];
	return(0);
}

/*****************************************************************
**	NCS Test: READING a file
*****************************************************************/

int test_openview( char *szInputFilename, BOOLEAN bRandomReads, BOOLEAN bReportTime )
{

	NCSFileView *pNCSFileView;
	NCSFileInfo	*pNCSFileInfo;

	NCSError eError = NCS_SUCCESS;
	UINT8	*p_output_buffer = NULL;
	UINT32	x_size, y_size, number_x, number_y;
	UINT32	start_x, start_y, end_x, end_y;
	int		regions;
	double	total_pixels = 0.0;
	UINT32	band;
	UINT32	nBands;
	UINT32	*band_list = NULL;	/* list of individual bands to read, may be subset of actual bands */
	INT32 nEPSG = -1;
	clock_t	start_time, mark_time;
	UINT32 nMaxWindow;

	printf("ECW READ EXAMPLE\n");
	start_time = mark_time = clock();
	/*
	**	Open the input NCSFileView
	*/
	eError = NCSOpenFileViewA(szInputFilename, &pNCSFileView, NULL);

	if (eError != NCS_SUCCESS) {
		printf("Could not open view for file:%s\n",szInputFilename);
		printf("Error = %s\n", NCSGetErrorText(eError));
		return(1);
	}
	NCSGetViewFileInfo(pNCSFileView, &pNCSFileInfo);
	x_size = pNCSFileInfo->nSizeX;
	y_size = pNCSFileInfo->nSizeY;
	nBands = pNCSFileInfo->nBands;
	printf("Input file is [%ld x %ld by %d bands]\n",
		(long)x_size, (long)y_size, nBands);
	nBands = 3;

	if(NCSGetEPSGCode(pNCSFileInfo->szProjection, pNCSFileInfo->szDatum, &nEPSG) == NCS_SUCCESS) {
		printf("Input file is georeferenced [EPSG:%d]\n", nEPSG);
	} else if(pNCSFileInfo->szProjection && pNCSFileInfo->szDatum) {
		printf("Input file is georeferenced [Projection:%s, Datum:%s]\n", pNCSFileInfo->szProjection, pNCSFileInfo->szDatum);
	} else {
		printf("Input file is not georeferenced\n");
	}
	fflush(stdout);

	// Have to set up the band list. Compatible with ER Mapper's method.
	// In this example we always request all bands.
	band_list = (UINT32 *) malloc(sizeof(UINT32) * nBands);
	if( !band_list ) {
		printf("Error - unable to malloc band list\n");
		NCSCloseFileView(pNCSFileView);
		return(1);
	}
	for( band = 0; band < nBands; band++ )
		band_list[band] = band;

// All of image
	start_x = 0;		start_y = 0;
	end_x = x_size - 1;	end_y = y_size - 1;
	number_x = x_size;	number_y = y_size;

	// error when MAX_WINDOW > image width or height
	nMaxWindow = MAX_WINDOW;
	nMaxWindow = NCSMin(nMaxWindow, x_size);
	nMaxWindow = NCSMin(nMaxWindow, y_size);

	{
		printf("[%d,%d] to [%d,%d] for [%d,%d]\n",
					 start_x, start_y, end_x, end_y, number_x, number_y);
		fflush(stdout);

		FILE *rgb_out = fopen("rgb.out", "w");

		eError = NCSSetFileView(pNCSFileView, 
						nBands, band_list,
						start_x, start_y, end_x, end_y,
						number_x, number_y);
		if( eError != NCS_SUCCESS) {
			printf("Error while setting file view to %d bands, TL[%d,%d] BR[%d,%d], Window size [%d,%d]\n",
				nBands, start_x, start_y, end_x, end_y, number_x, number_y);
			printf("Error = %s\n", NCSGetErrorText(eError));
			NCSCloseFileView(pNCSFileView);
			free(band_list);
			return(1);
		}

		UINT8 *rgb_line = (UINT8*)malloc(number_x * 3);
		int line;
		for (line = 0; line < number_y; line++) {
			NCSReadStatus eReadStatus;
			eReadStatus = NCSReadViewLineRGB( pNCSFileView, rgb_line);
			if (eReadStatus != NCS_READ_OK) {
				printf("Read line error at line %d\n",line);
				printf("Status code = %d\n", eReadStatus);
				NCSCloseFileView(pNCSFileView);
				return(1);
			}
			fwrite(rgb_line, 1, number_x * 3, rgb_out);
		}

		fflush(rgb_out);
		fclose(rgb_out);
		printf("lines: %d\n", line);
	}

	// Set bFreeCacheFile to TRUE if this is the last view and you want
	// to close the file, otherwise it will be kept open in the cache. 
	BOOLEAN bFreeCacheFile = FALSE;
	NCSCloseFileViewEx(pNCSFileView, bFreeCacheFile);
	free(band_list);

	return(0);
}
