#pragma comment(lib, "OpenCL.lib")
#include <CL/cl.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <setjmp.h>
#include <omp.h>
#include <time.h>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <io.h>

#define uchar unsigned char
#define uint unsigned int

extern "C"
{
#define JPEG_INTERNALS
#include "jpeglib.h"
#include "jdct.h"
//#include "jpegfdctfloat.c"
//#include "jpegidctfloat.c"
#include "jerror.h"
#include "jconfig.h"
#include "jmorecfg.h"
// void idct_simple(float *coef_block, float *output_buf, int stride);
// void fdct_simple(float *input_data_float, int stride, float * coef_block);
}

const int LOCAL_WINDOW_SIZE = 16;

jpeg_decompress_struct cinfo;
JSAMPROW input_image;
struct my_error_mgr {
	struct jpeg_error_mgr pub;
	jmp_buf setjmp_buffer;
};

typedef struct my_error_mgr * my_error_ptr;

METHODDEF(void)
my_error_exit(j_common_ptr cinfo)
{
	my_error_ptr myerr = (my_error_ptr)cinfo->err;
	(*cinfo->err->output_message) (cinfo);
	longjmp(myerr->setjmp_buffer, 1);
}

int read_img(const char *filename)
{
	struct my_error_mgr jerr;
	FILE *infile;
	JSAMPARRAY buffer;
	int row_stride;
	if ((infile = fopen(filename, "rb")) == NULL)
	{
		fprintf(stderr, "can't open %s\n", filename);
		return 0;
	}

	cinfo.err = jpeg_std_error(&jerr.pub);
	jerr.pub.error_exit = my_error_exit;

	if (setjmp(jerr.setjmp_buffer))
	{
		jpeg_destroy_decompress(&cinfo);
		fclose(infile);
		return 0;
	}

	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, infile);
	(void)jpeg_read_header(&cinfo, TRUE);
	cinfo.out_color_space = JCS_GRAYSCALE;
	(void)jpeg_start_decompress(&cinfo);
	row_stride = cinfo.output_width * cinfo.output_components;
	buffer = (*cinfo.mem->alloc_sarray)((j_common_ptr)&cinfo, JPOOL_IMAGE, row_stride, 1);
	input_image = new JSAMPLE[cinfo.image_height*cinfo.image_width];

	int index = 0;
	while (cinfo.output_scanline < cinfo.output_height)
	{
		(void)jpeg_read_scanlines(&cinfo, buffer, 1);
		for (size_t i = 0; i < row_stride; i += cinfo.output_components)
		{
			input_image[index++] = buffer[0][i];
		}
	}
	(void)jpeg_finish_decompress(&cinfo);
	fclose(infile);
	return 0;
}

uint getQuant(const uint num)
{
	return cinfo.quant_tbl_ptrs[0]->quantval[num];
}

void write_img(const char *filename, uchar *DCT_mtrx)
{
	FILE *file = fopen(filename, "wb");
	fprintf(file, "P5\n%i %i\n255\n", cinfo.image_width, cinfo.image_height);
	fwrite(matDCT, sizeof(uchar), cinfo.image_height * cinfo.image_width, file);
	fclose(file);
}

int main()
{
	cl_int err, status;
	cl_uint numDevices, numPlatforms;
	cl_device_id *devices, device;

	status = clGetPlatformIDs(0, NULL, &numPlatforms);
	//printf("Number of Platforms: %i\n", numPlatforms);

	cl_platform_id *platforms = new cl_platform_id[numPlatforms];
	status = clGetPlatformIDs(numPlatforms, platforms, NULL);

	for (int i = 0; i < (int)numPlatforms; i++)
	{
		size_t paramValSize(0);
		char *paramValue;

		clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 0, NULL, &paramValSize);
		paramValue = new char[paramValSize];
		clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, paramValSize, paramValue, NULL);
		//printf("Name of Platform: %s\n", str);
		delete[] paramValue;

		status = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
		//printf("\nNumber of Devices: %i\n", numDevices);
		devices = new cl_device_id[numDevices];
		status = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);

		for (int i = 0; i < (int)numDevices; i++)
		{
			size_t size;
			cl_device_type type;

			clGetDeviceInfo(devices[i], CL_DEVICE_TYPE, sizeof(type), &type, NULL);
			//printf("Type of Device N.%i: ", i + 1);

			switch (type)
			{
			case CL_DEVICE_TYPE_CPU:
				//printf("CPU(%d)\n", (int)type);
				break;
			case CL_DEVICE_TYPE_GPU:
				//rintf("GPU(%d)\n", (int)type);
				break;
			case CL_DEVICE_TYPE_ACCELERATOR:
				//printf("ACCELERATOR(%d)\n", (int)type);
				break;
			default:
				//printf("OTHER\n");
				break;
			}

			if (type == 2 || type == 4)
			{
				device = devices[i];
				goto found;
			}
		}
	}
	return -1;

found:

	cl_context context = clCreateContext(0, 1, &device, NULL, NULL, &status);
	cl_command_queue commandQueue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);

	FILE *fp;
	if (fp == NULL) {
		printf("File is not open\n");
		return 1;
	}
	fseek(fp, 0L, SEEK_END);
	size_t program_size = ftell(fp);
	fseek(fp, 0L, SEEK_SET);

	char* program_source = new char[program_size];
	fread(program_str, program_size, 1, fp);
	fclose(fp);

	program_size = _filelength(_fileno(fp));
	program_source = (char*)malloc(program_size + 1);
	fread(program_source, sizeof(char), program_size, fp);

	cl_program program = clCreateProgramWithSource(context, 1, (const char**)&program_source, &program_size, &status);
	status = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	if (program == NULL)
	{
		printf("\n[ERROR] Failed to read program.\n\n");
		return -1;
	}

	if (status != 0)
	{
		size_t len = 0;
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
		char *build = (char*)malloc(len);
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, len, build, NULL);
		printf("\n[ERROR] Kernel not compiled!\n");
		printf("\n\nBuild log:  %s\n\n", build);
		return -1;
	}

	cl_kernel mainKernel = clCreateKernel(program, "mainKernel", NULL);
	if (!mainKernel) {
		printf("\n[ERROR] Main kernel error.\n");
		return -1;
	}
	cl_kernel rightKernel = clCreateKernel(program, "rightKernel", NULL);
	if (!rightKernel) {
		printf("\n[ERROR] Right kernel error.\n");
		return -1;
	}
	cl_kernel downKernel = clCreateKernel(program, "downKernel", NULL);
	if (!downKernel) {
		printf("\n[ERROR] Right kernel error.\n");
		return -1;
	}
	cl_kernel cornerKernel = clCreateKernel(program, "cornerKernel", NULL);
	if (!angleKernel) {
		printf("[ERROR] Corner kernel error.\n");
		return -1;
	}
	cl_kernel division = clCreateKernel(program, "division", NULL);
	if (!division) {
		printf("[ERROR] Division kernel error.\n");
		return -1;
	}

	read_img("test.jpg");
	uint rows = cinfo.image_height;
	uint cols = cinfo.image_width;

	static const double aanscalefactor[DCTSIZE] = {
		1.0, 1.387039845, 1.306562965, 1.175875602,
		1.0, 0.785694958, 0.541196100, 0.275899379
	};

	float quant_matrix[DCTSIZE2];
	for (int i = 0; i<DCTSIZE; i++)
	{
		for (int j = 0; j<DCTSIZE; j++)
		{
			quant_matrix[i*DCTSIZE + j] = 1.0 / (aanscalefactor[i] * aanscalefactor[j] * 8.0 * getQuant(i*DCTSIZE + j));
		}
	}

	float dequant_matrix[DCTSIZE2];
	for (int i = 0; i<DCTSIZE; i++)
	{
		for (int j = 0; j<DCTSIZE; j++)
		{
			dequant_matrix[i*DCTSIZE + j] = (aanscalefactor[i] * aanscalefactor[j] * 0.125 * getQuant(i*DCTSIZE + j));
		}
	}

	cl_mem input_buff, DCT_buff, quant_mtrx_buff, dequant_mtrx_buff;
	input_buff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_uchar)*rows*cols, NULL, NULL);
	DCT_buff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_short)*rows*cols, NULL, NULL);
	dequant_mtrx_buff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float)*DCTSIZE2, NULL, NULL);
	quant_mtrx_buff = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float)*DCTSIZE2, NULL, NULL);

	if (!input_buff || !DCT_buff || !quant_mtrx_buff || !dequant_mtrx_buff)
	{
		printf("\n[ERROR] Memory allocation failed\n");
		return -1;
	}

	short *DCT_mtrx = new short[rows*cols];
	uint offset = 0;

	clEnqueueWriteBuffer(commandQueue, input_buff, CL_FALSE, 0, sizeof(cl_uchar)*rows*cols, input_image, 0, NULL, NULL);
	clEnqueueWriteBuffer(commandQueue, mem0, CL_FALSE, 0, sizeof(cl_short)*rows*cols, DCT_mtrx, 0, NULL, NULL);
	clEnqueueWriteBuffer(commandQueue, mem0, CL_FALSE, 0, sizeof(cl_float)*DCTSIZE2, quant_matrix, 0, NULL, NULL);
	clEnqueueWriteBuffer(commandQueue, mem0, CL_FALSE, 0, sizeof(cl_float)*DCTSIZE2, dequant_matrix, 0, NULL, NULL);

	clSetKernelArg(mainKernel, 0, sizeof(cl_mem), &input_buff);
	clSetKernelArg(mainKernel, 1, sizeof(cl_mem), &rows);
	clSetKernelArg(mainKernel, 2, sizeof(cl_mem), &cols);
	clSetKernelArg(mainKernel, 3, sizeof(cl_mem), &DCT_buff);
	clSetKernelArg(mainKernel, 4, sizeof(cl_mem), &quant_mtrx_buff);
	clSetKernelArg(mainKernel, 5, sizeof(cl_mem), &dequant_mtrx_buff);

	short tmp = 0;
	clEnqueueFillBuffer(commandQueue, DCT_buff, &tmp, 1, 0, sizeof(short)*rows*cols, 0, NULL, NULL);

	cl_event kernel_event;
	for (int y_offset = 0; y_offset <= DCTSIZE; y_offset += DCTSIZE)
	{
		for (int x_offset = 0; x_offset <= DCTSIZE; x_offset += DCTSIZE)
		{
			int y_start = y_offset;
			int y_end = ((rows - y_offset) / LOCAL_WINDOW_SIZE)*LOCAL_WINDOW_SIZE + y_offset;
			int x_start = x_offset;
			int x_end = ((cols - x_offset) / LOCAL_WINDOW_SIZE)*LOCAL_WINDOW_SIZE + x_offset;
			offset = y_offset*cols + x_offset;
			clSetKernelArg(mainKernel, 6, sizeof(cl_uint), &offset);
			size_t dct_local_work_size[2] = { DCTSIZE2, 1 };
			size_t newRows = ((y_end - y_start) / LOCAL_WINDOW_SIZE)*dct_local_work_size[1];
			size_t newCols = ((x_end - x_start) / LOCAL_WINDOW_SIZE)*dct_local_work_size[0];
			size_t dct_global_work_size[2] = { newCols, newRows };


			if (clEnqueueNDRangeKernel(command_queue, mainKernel, 2, NULL, dct_global_work_size, dct_local_work_size, 0, NULL, &kernel_event)) {
				printf("[ERROR] Main kernel execution failed.\n");
				return -1;
			}
		}
	}

	clSetKernelArg(rightKernel, 0, sizeof(cl_mem), &input_buff);
	clSetKernelArg(rightKernel, 1, sizeof(cl_mem), &rows);
	clSetKernelArg(rightKernel, 2, sizeof(cl_mem), &cols);
	clSetKernelArg(rightKernel, 3, sizeof(cl_mem), &DCT_buff);
	clSetKernelArg(rightKernel, 4, sizeof(cl_mem), &quant_mtrx_buff);
	clSetKernelArg(rightKernel, 5, sizeof(cl_mem), &dequant_mtrx_buff);

	size_t right_edge_DCT_local_size[1] = { DCTSIZE };
	size_t right_edge_DCT_global_size[1] = { rows / 2 };
	offset = 0;

	clSetKernelArg(rightKernel, 6, sizeof(cl_uint), &offset);

	if (clEnqueueNDRangeKernel(command_queue, rightKernel, 1, NULL, right_edge_DCT_global_size, ight_edge_DCT_local_size, 0, NULL, &kernel_event))
	{
		printf("\n[ERROR] Right kernel execution failed.\n");
		return -1;
	}

	right_edge_DCT_local_size[0] = DCTSIZE;
	right_edge_DCT_global_size[0] = (rows - LOCAL_WINDOW_SIZE) / 2;
	offset = cols*DCTSIZE;

	clSetKernelArg(rightKernel, 6, sizeof(cl_uint), &offset);

	if (clEnqueueNDRangeKernel(command_queue, rightKernel, 1, NULL, right_edge_DCT_global_size, ight_edge_DCT_local_size, 0, NULL, &kernel_event))
	{
		printf("\n[ERROR] Right kernel execution failed.\n");
		return -1;
	}

	clSetKernelArg(downKernel, 0, sizeof(cl_mem), &input_buff);
	clSetKernelArg(downKernel, 1, sizeof(cl_mem), &rows);
	clSetKernelArg(downKernel, 2, sizeof(cl_mem), &cols);
	clSetKernelArg(downKernel, 3, sizeof(cl_mem), &DCT_buff);
	clSetKernelArg(downKernel, 4, sizeof(cl_mem), &quant_mtrx_buff);
	clSetKernelArg(downKernel, 5, sizeof(cl_mem), &dequant_mtrx_buff);

	size_t down_edge_DCT_local_size[1] = { DCTSIZE };
	size_t down_edge_DCT_global_size[1] = { cols / 2 };
	offset = 0;

	clSetKernelArg(downKernel, 6, sizeof(cl_uint), &offset);

	if (clEnqueueNDRangeKernel(command_queue, downKernel, 1, NULL, down_edge_DCT_global_size, down_edge_DCT_local_size, 0, NULL, &kernel_event))
	{
		printf("\n[ERROR] Down kernel execution failed.\n");
		return -1;
	}

	down_edge_DCT_local_size[0] = DCTSIZE;
	down_edge_DCT_global_size[0] = (cols - LOCAL_WINDOW_SIZE) / 2;
	offset = DCTSIZE;

	clSetKernelArg(downKernel, 6, sizeof(cl_uint), &offset);

	if (clEnqueueNDRangeKernel(command_queue, downKernel, 1, NULL, down_edge_DCT_global_size, down_edge_DCT_local_size, 0, NULL, &kernel_event))
	{
		printf("\n[ERROR] Down kernel execution failed.\n");
		return -1;
	}

	clSetKernelArg(division, 0, sizeof(cl_mem), &input_buff);
	clSetKernelArg(division, 1, sizeof(cl_mem), &rows);
	clSetKernelArg(division, 2, sizeof(cl_mem), &cols);
	clSetKernelArg(division, 3, sizeof(cl_mem), &DCT_buff);
	clSetKernelArg(division, 4, sizeof(cl_mem), &quant_mtrx_buff);
	clSetKernelArg(division, 5, sizeof(cl_mem), &dequant_mtrx_buff);

	size_t global_size[2] = { cols / DCTSIZE, rows };

	if (clEnqueueNDRangeKernel(command_queue, division, 2, NULL, global_size, NULL, 0, NULL, &kernel_event))
	{
		printf("\n[ERROR] Division kernel execution failed.\n");
		return -1;
	}

	uchar *result = new uchar[rows * cols];
	clEnqueueReadBuffer(commandQueue, input_buff, TRUE, 0, sizeof(cl_uchar)*rows*cols, result, 0, NULL, NULL);

	write_img("result.pgm", result);
	delete[] DCT_mtrx;

	return 0;
}
