#define DCTSIZE 8
#define DCTSIZE2 64
#define LOCAL_BLOCK_SIZE 16
#define range_limit(x) (x)


#define uchar unsigned char
#define uint unsigned int

void jpeg_fdct_simple (float *output_data, float *input_data)
{
  float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  float tmp10, tmp11, tmp12, tmp13;
  float z1, z2, z3, z4, z5, z11, z13;
  int ctr;

  /* Pass 1: process rows. */

  for (ctr = 0; ctr < DCTSIZE; ctr++) {

    /* Load output_data into workspace */
  tmp0 = input_data[ctr*DCTSIZE + 0] + input_data[ctr*DCTSIZE + 7];
	tmp7 = input_data[ctr*DCTSIZE + 0] - input_data[ctr*DCTSIZE + 7];
	tmp1 = input_data[ctr*DCTSIZE + 1] + input_data[ctr*DCTSIZE + 6];
	tmp6 = input_data[ctr*DCTSIZE + 1] - input_data[ctr*DCTSIZE + 6];
	tmp2 = input_data[ctr*DCTSIZE + 2] + input_data[ctr*DCTSIZE + 5];
	tmp5 = input_data[ctr*DCTSIZE + 2] - input_data[ctr*DCTSIZE + 5];
	tmp3 = input_data[ctr*DCTSIZE + 3] + input_data[ctr*DCTSIZE + 4];
	tmp4 = input_data[ctr*DCTSIZE + 3] - input_data[ctr*DCTSIZE + 4];

    /* Even part */

    tmp10 = tmp0 + tmp3;  /* phase 2 */
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;

    /* Apply unsigned->signed conversion */
    output_data[ctr*DCTSIZE + 0] = tmp10 + tmp11 - 8 * 128; /* phase 3 */
    output_data[ctr*DCTSIZE + 4] = tmp10 - tmp11;

    z1 = (tmp12 + tmp13) * 0.707106781f; /* c4 */
    output_data[ctr*DCTSIZE + 2] = tmp13 + z1;  /* phase 5 */
    output_data[ctr*DCTSIZE + 6] = tmp13 - z1;

    /* Odd part */

    tmp10 = tmp4 + tmp5;  /* phase 2 */
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;

    /* The rotator is modified from fig 4-8 to avoid extra negations. */
    z5 = (tmp10 - tmp12) * 0.382683433f; /* c6 */
    z2 = 0.541196100f * tmp10 + z5; /* c2-c6 */
    z4 = 1.306562965f * tmp12 + z5; /* c2+c6 */
    z3 = tmp11 * 0.707106781f; /* c4 */

    z11 = tmp7 + z3;    /* phase 5 */
    z13 = tmp7 - z3;

    output_data[ctr*DCTSIZE + 5] = z13 + z2;  /* phase 6 */
    output_data[ctr*DCTSIZE + 3] = z13 - z2;
    output_data[ctr*DCTSIZE + 1] = z11 + z4;
    output_data[ctr*DCTSIZE + 7] = z11 - z4;

       /* advance pointer to next row */
  }

  /* Pass 2: process columns. */

  for (ctr = 0; ctr < DCTSIZE; ctr++) {

    tmp0 = output_data[DCTSIZE * 0 + ctr] + output_data[DCTSIZE * 7 + ctr];
    tmp7 = output_data[DCTSIZE * 0 + ctr] - output_data[DCTSIZE * 7 + ctr];
    tmp1 = output_data[DCTSIZE * 1 + ctr] + output_data[DCTSIZE * 6 + ctr];
    tmp6 = output_data[DCTSIZE * 1 + ctr] - output_data[DCTSIZE * 6 + ctr];
    tmp2 = output_data[DCTSIZE * 2 + ctr] + output_data[DCTSIZE * 5 + ctr];
    tmp5 = output_data[DCTSIZE * 2 + ctr] - output_data[DCTSIZE * 5 + ctr];
    tmp3 = output_data[DCTSIZE * 3 + ctr] + output_data[DCTSIZE * 4 + ctr];
    tmp4 = output_data[DCTSIZE * 3 + ctr] - output_data[DCTSIZE * 4 + ctr];

    /* Even part */

    tmp10 = tmp0 + tmp3;  /* phase 2 */
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;

    output_data[DCTSIZE * 0 + ctr] = tmp10 + tmp11; /* phase 3 */
    output_data[DCTSIZE * 4 + ctr] = tmp10 - tmp11;

    z1 = (tmp12 + tmp13) * 0.707106781f; /* c4 */
    output_data[DCTSIZE * 2 + ctr] = tmp13 + z1; /* phase 5 */
    output_data[DCTSIZE * 6 + ctr] = tmp13 - z1;

    /* Odd part */

    tmp10 = tmp4 + tmp5;  /* phase 2 */
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;

    /* The rotator is modified from fig 4-8 to avoid extra negations. */
    z5 = (tmp10 - tmp12) * 0.382683433f; /* c6 */
    z2 = 0.541196100f * tmp10 + z5; /* c2-c6 */
    z4 = 1.306562965f * tmp12 + z5; /* c2+c6 */
    z3 = tmp11 * 0.707106781f; /* c4 */

    z11 = tmp7 + z3;    /* phase 5 */
    z13 = tmp7 - z3;

    output_data[DCTSIZE * 5 + ctr] = z13 + z2; /* phase 6 */
    output_data[DCTSIZE * 3 + ctr] = z13 - z2;
    output_data[DCTSIZE * 1 + ctr] = z11 + z4;
    output_data[DCTSIZE * 7 + ctr] = z11 - z4;

  }
}

void jpeg_idct_simple(const float constant *quantptr, float *coef_block)
{
  float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  float tmp10, tmp11, tmp12, tmp13;
  float z5, z10, z11, z12, z13;

  /* Pass 1: process columns from input, store into work array. */

  for (int ctr = 0; ctr < DCTSIZE; ctr++) {


    /* Even part */

    tmp0 = coef_block[DCTSIZE * 0 + ctr] * quantptr[DCTSIZE * 0 + ctr];
    tmp1 = coef_block[DCTSIZE * 2 + ctr] * quantptr[DCTSIZE * 2 + ctr];
    tmp2 = coef_block[DCTSIZE * 4 + ctr] * quantptr[DCTSIZE * 4 + ctr];
    tmp3 = coef_block[DCTSIZE * 6 + ctr] * quantptr[DCTSIZE * 6 + ctr];



    tmp10 = tmp0 + tmp2;  /* phase 3 */
    tmp11 = tmp0 - tmp2;

    tmp13 = tmp1 + tmp3;  /* phases 5-3 */
    tmp12 = (tmp1 - tmp3) * 1.414213562f - tmp13; /* 2*c4 */

    tmp0 = tmp10 + tmp13;  /* phase 2 */
    tmp3 = tmp10 - tmp13;
    tmp1 = tmp11 + tmp12;
    tmp2 = tmp11 - tmp12;

    /* Odd part */

	tmp4 = coef_block[DCTSIZE * 1 + ctr] * quantptr[DCTSIZE * 1 + ctr];
    tmp5 = coef_block[DCTSIZE * 3 + ctr] * quantptr[DCTSIZE * 3 + ctr];
    tmp6 = coef_block[DCTSIZE * 5 + ctr] * quantptr[DCTSIZE * 5 + ctr];
    tmp7 = coef_block[DCTSIZE * 7 + ctr] * quantptr[DCTSIZE * 7 + ctr];


    z13 = tmp6 + tmp5;    /* phase 6 */
    z10 = tmp6 - tmp5;
    z11 = tmp4 + tmp7;
    z12 = tmp4 - tmp7;

    tmp7 = z11 + z13;    /* phase 5 */
    tmp11 = (z11 - z13) * 1.414213562f; /* 2*c4 */

    z5 = (z10 + z12) * 1.847759065f; /* 2*c2 */
    tmp10 = z5 - z12 * 1.082392200f; /* 2*(c2-c6) */
    tmp12 = z5 - z10 * 2.613125930f; /* 2*(c2+c6) */

    tmp6 = tmp12 - tmp7;  /* phase 2 */
    tmp5 = tmp11 - tmp6;
    tmp4 = tmp10 - tmp5;

    coef_block[DCTSIZE * 0 + ctr] = tmp0 + tmp7;
    coef_block[DCTSIZE * 7 + ctr] = tmp0 - tmp7;
    coef_block[DCTSIZE * 1 + ctr] = tmp1 + tmp6;
    coef_block[DCTSIZE * 6 + ctr] = tmp1 - tmp6;
    coef_block[DCTSIZE * 2 + ctr] = tmp2 + tmp5;
    coef_block[DCTSIZE * 5 + ctr] = tmp2 - tmp5;
    coef_block[DCTSIZE * 3 + ctr] = tmp3 + tmp4;
    coef_block[DCTSIZE * 4 + ctr] = tmp3 - tmp4;

  }

  /* Pass 2: process rows from work array, store into output array. */


  for (int ctr = 0; ctr < DCTSIZE; ctr++) {

    /* Even part */

    /* Apply signed->unsigned and prepare float->int conversion */
    z5 = coef_block[ctr*DCTSIZE + 0] + (128 + 0.5f);
    tmp10 = z5 + coef_block[ctr*DCTSIZE + 4];
    tmp11 = z5 - coef_block[ctr*DCTSIZE + 4];

    tmp13 = coef_block[ctr*DCTSIZE + 2] + coef_block[ctr*DCTSIZE + 6];
    tmp12 = (coef_block[ctr*DCTSIZE + 2] - coef_block[ctr*DCTSIZE + 6]) * 1.414213562f - tmp13;

    tmp0 = tmp10 + tmp13;
    tmp3 = tmp10 - tmp13;
    tmp1 = tmp11 + tmp12;
    tmp2 = tmp11 - tmp12;

    /* Odd part */

    z13 = coef_block[ctr*DCTSIZE + 5] + coef_block[ctr*DCTSIZE + 3];
    z10 = coef_block[ctr*DCTSIZE + 5] - coef_block[ctr*DCTSIZE + 3];
    z11 = coef_block[ctr*DCTSIZE + 1] + coef_block[ctr*DCTSIZE + 7];
    z12 = coef_block[ctr*DCTSIZE + 1] - coef_block[ctr*DCTSIZE + 7];

    tmp7 = z11 + z13;
    tmp11 = (z11 - z13) * 1.414213562f;

    z5 = (z10 + z12) * 1.847759065f; /* 2*c2 */
    tmp10 = z5 - z12 * 1.082392200f; /* 2*(c2-c6) */
    tmp12 = z5 - z10 * 2.613125930f; /* 2*(c2+c6) */

    tmp6 = tmp12 - tmp7;
    tmp5 = tmp11 - tmp6;
    tmp4 = tmp10 - tmp5;

    /* Final output stage: float->int conversion and range-limit */

    coef_block[ctr*DCTSIZE + 0] = range_limit(tmp0 + tmp7);
    coef_block[ctr*DCTSIZE + 7] = range_limit(tmp0 - tmp7);
    coef_block[ctr*DCTSIZE + 1] = range_limit(tmp1 + tmp6);
    coef_block[ctr*DCTSIZE + 6] = range_limit(tmp1 - tmp6);
    coef_block[ctr*DCTSIZE + 2] = range_limit(tmp2 + tmp5);
    coef_block[ctr*DCTSIZE + 5] = range_limit(tmp2 - tmp5);
    coef_block[ctr*DCTSIZE + 3] = range_limit(tmp3 + tmp4);
    coef_block[ctr*DCTSIZE + 4] = range_limit(tmp3 - tmp4);
  }
}

void quantization(float *input_data, float *output_data,  constant const float  *quant_matrix,  constant const float  *dequant_matrix)
{

  jpeg_fdct_simple(output_data, input_data);

  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      output_data[i * 8 + j] = rint(output_data[i * 8 + j] * quant_matrix[i * 8 + j]);
    }
  }

  jpeg_idct_simple(dequant_matrix, output_data);
}


kernel void mainKernel(const uchar global *input_img, const uint rows, const uint cols, short global *resDCT,  constant const float  *quant_matrix,  constant const float  *dequant_matrix, const uint offset)
{

  local float local_matrix[256];

  int x = get_local_id(0) % 8;
  int y = get_local_id(0) / 8;

  for (int i = 0; i <= 8; i += 8)
  {
    for (int j = 0; j <= 8; j += 8)
    {
      local_matrix[(y + i) * 16 + x + j] = convert_float(*(global uchar*)&input_img[(get_group_id(1) * LOCAL_BLOCK_SIZE + y + i) * cols + get_group_id(0) * LOCAL_BLOCK_SIZE + x + j + offset]);
    }
  }

  barrier(CLK_LOCAL_MEM_FENCE);

  float res[DCTSIZE2];
  float block8x8[DCTSIZE2];
  for (int i = 0; i < DCTSIZE; i++)
  {
    for (int j = 0; j < DCTSIZE; j++)
    {
      block8x8[i * DCTSIZE+j] = local_matrix[(y + i) * LOCAL_BLOCK_SIZE + x + j];
    }
  }
  quantization(block8x8, res, quant_matrix, dequant_matrix);

  barrier(CLK_LOCAL_MEM_FENCE);

  for (int i = 0; i <= 8; i += 8)
  {
    for (int j = 0; j <= 8; j += 8)
    {
      local_matrix[(y + i) * 16 + x + j] = 0.0f;
    }
  }

  barrier(CLK_LOCAL_MEM_FENCE);

  for (int i = 0; i < DCTSIZE; i++)
  {
    for (int j = 0; j < DCTSIZE; j++)
    {
      local_matrix[(y + i) * 16 + x + j] += res[i * DCTSIZE + j];
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }

  for (int i = 0; i <= 8; i += 8)
  {
    for (int j = 0; j <= 8; j += 8)
    {
      *(global short*)&resDCT[(get_group_id(1) * LOCAL_BLOCK_SIZE + y + i) * cols + get_group_id(0) * LOCAL_BLOCK_SIZE + x + j + offset] += convert_short_sat_rte(local_matrix[(y + i) * 16 + x + j]);
    }
  }

}


kernel void rightKernel(const uchar global *input_img, const uint rows, const uint cols, short global *resDCT,  constant const float  *quant_matrix,  constant const float  *dequant_matrix, const uint offset)
{
  local float local_matrix[LOCAL_BLOCK_SIZE*DCTSIZE];
  int x = get_local_id(0);

  for (int i = 0; i <= DCTSIZE; i += DCTSIZE)
  {
    for (int j = 0; j < DCTSIZE; j++)
    {
      local_matrix[(x + i) * DCTSIZE + j] = convert_float(*(global uchar*)&input_img[(get_group_id(0) * LOCAL_BLOCK_SIZE + x + i) * cols + j + cols - DCTSIZE + offset]);
    }
  }

  barrier(CLK_LOCAL_MEM_FENCE);

  float res[DCTSIZE2];
  float block8x8[DCTSIZE2];
  for (int i = 0; i < DCTSIZE; i++)
  {
    for (int j = 0; j < DCTSIZE; j++)
    {
      block8x8[i * DCTSIZE + j] = local_matrix[(x + i) * DCTSIZE + j];
    }
  }
  quantization(block8x8, res, quant_matrix, dequant_matrix);
  barrier(CLK_LOCAL_MEM_FENCE);


  for (int i = 0; i <= DCTSIZE; i += DCTSIZE)
  {
    for (int j = 0; j < DCTSIZE; j++)
    {
      local_matrix[(x + i) * DCTSIZE + j] = 0;
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  for (int i = 0; i < DCTSIZE; i++)
  {
    for (int j = 0; j < DCTSIZE; j++)
    {
      local_matrix[(x + i) * DCTSIZE + j] += res[i * DCTSIZE + j];
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }

  for (int i = 0; i <= DCTSIZE; i += DCTSIZE)
  {
    for (int j = 0; j < DCTSIZE; j++)
    {
      *(global short*)&resDCT[(get_group_id(0) * LOCAL_BLOCK_SIZE + x + i) * cols + j + cols - DCTSIZE + offset] += convert_short_sat_rte(local_matrix[(x + i) * DCTSIZE + j]);
    }
  }
}



kernel void downKernel(const uchar global *input_img, const uint rows, const uint cols, short global *resDCT,  constant const float  *quant_matrix,  constant const float  *dequant_matrix, const uint offset)
{
  local float local_matrix[DCTSIZE * LOCAL_BLOCK_SIZE];
  int x = get_local_id(0);

  for (int i = 0; i < LOCAL_BLOCK_SIZE; i++)
  {
    local_matrix[x * LOCAL_BLOCK_SIZE+i] = convert_float(*(global uchar*)&input_img[(rows - DCTSIZE+x) * cols + i + (get_group_id(0) * LOCAL_BLOCK_SIZE) + offset]);
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  float res[DCTSIZE2];
  float block8x8[DCTSIZE2];
  for (int i = 0; i < DCTSIZE; i++)
  {
    for (int j = 0; j < DCTSIZE; j++)
    {
      block8x8[i * DCTSIZE+j] = local_matrix[(x + i) * LOCAL_BLOCK_SIZE + j];
    }
  }
  quantization(block8x8, res, quant_matrix, dequant_matrix);
  barrier(CLK_LOCAL_MEM_FENCE);

  for (int i = 0; i<LOCAL_BLOCK_SIZE; i++)
  {
    local_matrix[x * LOCAL_BLOCK_SIZE + i] = 0;
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  for (int i = 0; i < DCTSIZE; i++)
  {
    for (int j = 0; j < DCTSIZE; j++)
    {
      local_matrix[i * LOCAL_BLOCK_SIZE + x + j] += res[i * DCTSIZE + j];
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }

  for (int i = 0; i < LOCAL_BLOCK_SIZE; i++)
  {
    *(global short*)&resDCT[(rows - DCTSIZE + x) * cols + i + (get_group_id(0) * LOCAL_BLOCK_SIZE + offset)] += convert_short_sat_rte(local_matrix[x * LOCAL_BLOCK_SIZE + i]);
  }

}

kernel void cornerKernel(const uchar global *input_img, const uint rows, const uint cols, short global *resDCT,  constant const float  *quant_matrix,  constant const float  *dequant_matrix)
{
  local float local_matrix[DCTSIZE2];
  int x = get_local_id(0);

  for (int i = 0; i < DCTSIZE; i++)
  {
    local_matrix[x * DCTSIZE + i] = convert_float(*(global uchar*)&input_img[(rows - DCTSIZE+x) * cols + i + cols - DCTSIZE]);
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  float res[DCTSIZE2];
  if (x == 0)
  {
    float block8x8[DCTSIZE2];
    for (int i = 0; i < DCTSIZE; i++)
    {
      for (int j = 0; j < DCTSIZE; j++)
      {
        block8x8[i * DCTSIZE + j] = local_matrix[(x + i) * DCTSIZE + j];
      }
    }
    quantization(block8x8, res, quant_matrix, dequant_matrix);
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  for (int i = 0; i < DCTSIZE; i++)
  {
    local_matrix[x * DCTSIZE + i] = 0;
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  if (x == 0) {
    for (int i = 0; i < DCTSIZE2; i++){
      local_matrix[i] = res[i];
    }
  }

  for (int i = 0; i < DCTSIZE; i++)
  {
    *(global short*)&resDCT[(rows - DCTSIZE + x) * cols + i + cols - DCTSIZE] += convert_short_sat_rte(local_matrix[x * DCTSIZE + i]);
  }
}

constant const float8 MatUL[DCTSIZE] = {
  (float8)(1/1.0f,  1/2.0f,  1/3.0f,  1/4.0f,  1/5.0f,  1/6.0f,  1/7.0f,  1/8.0f),
  (float8)(1/2.0f,  1/4.0f,  1/6.0f,  1/8.0f,  1/10.0f,  1/12.0f,  1/14.0f,  1/16.0f),
  (float8)(1/3.0f,  1/6.0f,  1/9.0f,  1/12.0f,  1/15.0f,  1/18.0f,  1/21.0f,  1/24.0f),
  (float8)(1/4.0f,  1/8.0f,  1/12.0f,  1/16.0f,  1/20.0f,  1/24.0f,  1/28.0f,  1/32.0f),
  (float8)(1/5.0f,  1/10.0f,  1/15.0f,  1/20.0f,  1/25.0f,  1/30.0f,  1/35.0f,  1/40.0f),
  (float8)(1/6.0f,  1/12.0f,  1/18.0f,  1/24.0f,  1/30.0f,  1/36.0f,  1/42.0f,  1/48.0f),
  (float8)(1/7.0f,  1/14.0f,  1/21.0f,  1/28.0f,  1/35.0f,  1/42.0f,  1/49.0f,  1/56.0f),
  (float8)(1/8.0f,  1/16.0f,  1/24.0f,  1/32.0f,  1/40.0f,  1/48.0f,  1/56.0f,  1/64.0f)
};

constant const float8 MatDL[DCTSIZE] = {
  (float8)(1/8.0f,  1/16.0f,  1/24.0f,  1/32.0f,  1/40.0f,  1/48.0f,  1/56.0f,  1/64.0f),
  (float8)(1/7.0f,  1/14.0f,  1/21.0f,  1/28.0f,  1/35.0f,  1/42.0f,  1/49.0f,  1/56.0f),
  (float8)(1/6.0f,  1/12.0f,  1/18.0f,  1/24.0f,  1/30.0f,  1/36.0f,  1/42.0f,  1/48.0f),
  (float8)(1/5.0f,  1/10.0f,  1/15.0f,  1/20.0f,  1/25.0f,  1/30.0f,  1/35.0f,  1/40.0f),
  (float8)(1/4.0f,  1/8.0f,  1/12.0f,  1/16.0f,  1/20.0f,  1/24.0f,  1/28.0f,  1/32.0f),
  (float8)(1/3.0f,  1/6.0f,  1/9.0f,  1/12.0f,  1/15.0f,  1/18.0f,  1/21.0f,  1/24.0f),
  (float8)(1/2.0f,  1/4.0f,  1/6.0f,  1/8.0f,  1/10.0f,  1/12.0f,  1/14.0f,  1/16.0f),
  (float8)(1/1.0f,  1/2.0f,  1/3.0f,  1/4.0f,  1/5.0f,  1/6.0f,  1/7.0f,  1/8.0f)
};

constant const float8 MatUR[DCTSIZE] = {
  (float8)(1/8.0f,  1/7.0f,  1/6.0f,  1/5.0f,  1/4.0f,  1/3.0f,  1/2.0f,  1/1.0f),
  (float8)(1/16.0f,  1/14.0f,  1/12.0f,  1/10.0f,  1/8.0f,  1/6.0f,  1/4.0f,  1/2.0f),
  (float8)(1/24.0f,  1/21.0f,  1/18.0f,  1/15.0f,  1/12.0f,  1/9.0f,  1/6.0f,  1/3.0f),
  (float8)(1/32.0f,  1/28.0f,  1/24.0f,  1/20.0f,  1/16.0f,  1/12.0f,  1/8.0f,  1/4.0f),
  (float8)(1/40.0f,  1/35.0f,  1/30.0f,  1/25.0f,  1/20.0f,  1/15.0f,  1/10.0f,  1/5.0f),
  (float8)(1/48.0f,  1/42.0f,  1/36.0f,  1/30.0f,  1/24.0f,  1/18.0f,  1/12.0f,  1/6.0f),
  (float8)(1/56.0f,  1/49.0f,  1/42.0f,  1/35.0f,  1/28.0f,  1/21.0f,  1/14.0f,  1/7.0f),
  (float8)(1/64.0f,  1/56.0f,  1/48.0f,  1/40.0f,  1/32.0f,  1/24.0f,  1/16.0f,  1/8.0f)
};

constant const float8 MatDR[DCTSIZE] = {
  (float8)(1/64.0f,  1/56.0f,  1/48.0f,  1/40.0f,  1/32.0f,  1/24.0f,  1/16.0f,  1/8.0f),
  (float8)(1/56.0f,  1/49.0f,  1/42.0f,  1/35.0f,  1/28.0f,  1/21.0f,  1/14.0f,  1/7.0f),
  (float8)(1/48.0f,  1/42.0f,  1/36.0f,  1/30.0f,  1/24.0f,  1/18.0f,  1/12.0f,  1/6.0f),
  (float8)(1/40.0f,  1/35.0f,  1/30.0f,  1/25.0f,  1/20.0f,  1/15.0f,  1/10.0f,  1/5.0f),
  (float8)(1/32.0f,  1/28.0f,  1/24.0f,  1/20.0f,  1/16.0f,  1/12.0f,  1/8.0f,  1/4.0f),
  (float8)(1/24.0f,  1/21.0f,  1/18.0f,  1/15.0f,  1/12.0f,  1/9.0f,  1/6.0f,  1/3.0f),
  (float8)(1/16.0f,  1/14.0f,  1/12.0f,  1/10.0f,  1/8.0f,  1/6.0f,  1/4.0f,  1/2.0f),
  (float8)(1/8.0f,  1/7.0f,  1/6.0f,  1/5.0f,  1/4.0f,  1/3.0f,  1/2.0f,  1/1.0f)
};

kernel void division(short global *resDCT, uchar global *result_img, const uint rows, const uint cols)
{
  int x = get_global_id(0);
  int y = get_global_id(1);

  float8 tmp = convert_float8(*(global short8*)&resDCT[y * cols + x * DCTSIZE]);

  if (y >= DCTSIZE && (x * DCTSIZE) >= DCTSIZE && y < (rows - DCTSIZE) && (x * DCTSIZE) < (cols - DCTSIZE)) {
    tmp *= (float8)(1/64.0f,  1/64.0f,  1/64.0f,  1/64.0f,  1/64.0f,  1/64.0f,  1/64.0f,  1/64.0f);
  } else if (y < DCTSIZE && (x * DCTSIZE) < DCTSIZE)
  {
    tmp *= MatUL[y % DCTSIZE];
  } else if (y >= (rows - DCTSIZE) && (x * DCTSIZE) < DCTSIZE)
  {
    tmp *= MatDL[y % DCTSIZE];
  } else if (y >= (rows - DCTSIZE) && (x * DCTSIZE) >= (cols - DCTSIZE))
  {
    tmp *= MatDR[y%DCTSIZE];
  } else if (y < DCTSIZE && (x * DCTSIZE) >= (cols - DCTSIZE))
  {
    tmp *= MatUR [y % DCTSIZE];
  } else if (y < DCTSIZE && (x * DCTSIZE) >= DCTSIZE && (x * DCTSIZE) < (cols - DCTSIZE))
  {
    tmp *= (float8)(1.0f / (((y % DCTSIZE) + 1) * DCTSIZE));
  } else if (y >= (rows - DCTSIZE) && (x * DCTSIZE) >= DCTSIZE && (x * DCTSIZE) < (cols - DCTSIZE))
  {
    tmp *= (float8)(1.0f / ((DCTSIZE - (y % DCTSIZE)) * DCTSIZE));
  } else if ((x * DCTSIZE) >= (cols - DCTSIZE) && y >= (DCTSIZE) && y < (rows - DCTSIZE))
  {
    tmp *= MatUR[DCTSIZE - 1];
  } else if ((x * DCTSIZE) < DCTSIZE && y >= (DCTSIZE) && y < (rows - DCTSIZE))
  {
    tmp *= MatUL[DCTSIZE - 1];
  }

  *(global uchar8*)&result_img[y * cols + x * DCTSIZE] = convert_uchar8_sat_rte(tmp);
}
