
--------------------------------------------------------------------------------------------------


Image Melding 
MATLAB Scripts and Mex Executables Version 1.0 (2012-03-12)

  By Soheil Darabi
  Copyright:
  
	This software is distributed for free for non-commercial use only.

	Copyright © 2012 Adobe Systems Inc.

	This software is provided by the copyright holders and the contributors 
	``as is'' and any express or implied warranties, including, but not limited 
	to, the implied warranties of merchantability and fitness for a particular 
	purpose are disclaimed. In no event shall the copyright holders or 
	contributors be liable for any direct, indirect, incidental, special, 
	exemplary, or consequential damages (including, but not limited to, 
	procurement of substitute goods or services; loss of use, data, or profits;
	or business interruption) however caused and on any theory of liability, 
	whether in contract, strict liability, or tort (including negligence or 
	otherwise) arising in any way out of the use of this software, even if 
	advised of the possibility of such damage.

--------------------------------------------------------------------------------------------------
Background
--------------------------------------------------------------------------------------------------

For information please see the paper:

 - Image Melding: Combining Inconsistent Images using Patch-based Synthesis
   SIGGRAPH 2012, Soheil Darabi, Eli Shechtman, Connelly Barnes, Dan B Goldman, Pradeep Sen
   http://agl.unm.edu/melding



Please cite this paper if you use this code in an academic publication.

--------------------------------------------------------------------------------------------------
Contents
--------------------------------------------------------------------------------------------------

 - We include mex files of search and vote part of the Image Melding algorithm in "Mexfiles"
	directory.  These mex files are un-optimized, slower than the implementation used for 
	the paper result, but are general and can be tunned to generate the best qulaity.
 - To understand how to use the mex file, we have also included MATLAB scripts for two 
	applications we demonstrated in the paper:
		1) HoleFilling.m: that shows how image melding can be used for the hole-filling application.
		2) TextureInterpolation.m:  that shows texture interpolation application of image melding.
		3) ImageCloning.m: that shows the cloning application of the paper.  The inputs are source,
			target, clone mask, and synthesis mask.  The algorithm copies part of the source , 
			according to clone mask, onto target.  The algorithm takes an additional input, the synthesis
			mask, which will be synthesized to generate a better blend between two sources.
		
		Other files:
		1) do_em_* these files are used to do the core synthesis iterations of the algorithm 
			and according to the application they are namded differently.
		2) do_poisson_iterations: solve the simple screen Poisson equation.  We just released 
			the slower linear equation solver using Matlab source codes.  To further accelerate
			the algorithm, you can use "Convol Pyarmids" source code released by Farbman et al.
			at: http://www.cs.huji.ac.il/labs/cglab/projects/convpyr/

If you find any bugs or have comments/questions, please contact 
Soheil Darabi at sdarabi@gmail.com.

--------------------------------------------------------------------------------------------------
How to use it
--------------------------------------------------------------------------------------------------

Depending to your application choose either of HoleFilling, TextureInterpolation,
 or ImageCloning and then find the name of the input image and assign it to your 
 own file and run the Matlab code.  The intermediate images produced by Melding 
 iteration will be stored in Results folder and the final result would be the last
 image in that folder.