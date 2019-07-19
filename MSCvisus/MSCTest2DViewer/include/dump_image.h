#include <cstdlib>


using namespace std;



//void write_image(string& filename, int wx, int wy) {
//   cout << wx << " " << wy << endl;
//
//   FILE* outfile = fopen(filename.c_str(), "wb");
//
//   if (!outfile) {
//      cerr << "Failed to open: " << filename << endl;
//      exit(1);
//   }
//
//
//   JSAMPLE* buffer = new JSAMPLE[3*wx*wy];
//   glReadPixels( 0,0,wx,wy, GL_RGB, GL_UNSIGNED_BYTE, buffer);
//
//   //mirror y axis
//   /*
//   for (int i=0; i<wx; i++) {
//      for (int j=0; j<wy; j++) {
//         int p0 = 3*(i*wy + j);
//         int p1 = 3*(i*wy + (wy - j));
//
//         for (int k = 0; k<3; k++) {
//            JSAMPLE c = buffer[p0];
//            buffer[p0] = buffer[p1];
//            buffer[p1] = c;
//            p0++;
//            p1++;
//         }
//      }
//   }
//   */
//
//
//
//
//   struct jpeg_compress_struct cinfo;
//   struct jpeg_error_mgr       jerr;
//
//   cinfo.err = jpeg_std_error(&jerr);
//   jpeg_create_compress(&cinfo);
//   jpeg_stdio_dest(&cinfo, outfile);
//
//   cinfo.image_width      = wx;
//   cinfo.image_height     = wy;
//   cinfo.input_components = 3;
//   cinfo.in_color_space   = JCS_RGB;
//
//   jpeg_set_defaults(&cinfo);
//   jpeg_set_quality (&cinfo, 100, true);
//   jpeg_start_compress(&cinfo, true);
//
//
//   JSAMPROW row_pointer;          
//   for (int j=0; j<wy; j++) {
//      row_pointer = (JSAMPROW) &buffer[(wy-cinfo.next_scanline-1)*3*wx];
//      jpeg_write_scanlines(&cinfo, &row_pointer, 1);
//   }
//
//   jpeg_finish_compress(&cinfo);
//
//   fclose(outfile);
// 
//   delete[] buffer;
//
//}
