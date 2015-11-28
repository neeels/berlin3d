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
int ecw2rgb(const char *in_path, const char *out_path);

int main(int argc, char **argv)
{

	char	*in_path;
	char	*out_path;
	int		nError = 0;

	/*
	 * 	Initialize the library if we are linking statically
	 */
	NCSInit();
	
	if (argc != 3) {
	  printf("Usage: %s <in.ecw> <out.rgb>\n", argv[0]);
	  return(1);
	}

	in_path = argv[1];
	out_path = argv[2];


	nError = ecw2rgb(in_path, out_path);
	if( nError ) {
		printf("Openview test returned an error\n");
		NCSShutdown();        
		return(nError);
	}

	NCSShutdown();

	return(0);
}



/*****************************************************************
**	NCS Test: READING a file
*****************************************************************/

int ecw2rgb(const char *in_path, const char *out_path)
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

	printf("ECW READ EXAMPLE\n");
	start_time = mark_time = clock();
	/*
	**	Open the input NCSFileView
	*/
	eError = NCSOpenFileViewA(in_path, &pNCSFileView, NULL);

	if (eError != NCS_SUCCESS) {
		printf("Could not open view for file:%s\n",in_path);
		printf("Error = %s\n", NCSGetErrorText(eError));
		return(1);
	}
	NCSGetViewFileInfo(pNCSFileView, &pNCSFileInfo);
	x_size = pNCSFileInfo->nSizeX;
	y_size = pNCSFileInfo->nSizeY;
	//nBands = pNCSFileInfo->nBands;
	nBands = 3;

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

	number_x = number_y = 8192;
	//printf("SCALING to %dx%d\n", number_x, number_y);

	{
		printf("[%d,%d] to [%d,%d] for [%d,%d]\n",
					 start_x, start_y, end_x, end_y, number_x, number_y);
		fflush(stdout);

		FILE *rgb_out;
		if (out_path[0] == '-' && out_path[1] == '\0')
			rgb_out = stdout;
		else
			rgb_out = fopen(out_path, "w");

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
		if (rgb_out != stdout)
			fclose(rgb_out);
	}

	// Set bFreeCacheFile to TRUE if this is the last view and you want
	// to close the file, otherwise it will be kept open in the cache. 
	BOOLEAN bFreeCacheFile = FALSE;
	NCSCloseFileViewEx(pNCSFileView, bFreeCacheFile);
	free(band_list);

	return(0);
}
