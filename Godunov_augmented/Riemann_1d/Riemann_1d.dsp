# Microsoft Developer Studio Project File - Name="Riemann_1d" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=Riemann_1d - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Riemann_1d.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Riemann_1d.mak" CFG="Riemann_1d - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Riemann_1d - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "Riemann_1d - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "Riemann_1d - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "Riemann_1d - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /imsl /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "Riemann_1d - Win32 Release"
# Name "Riemann_1d - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\barotropic.f90
DEP_F90_BAROT=\
	".\Debug\conserved_quantities.mod"\
	".\Debug\my_messages.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\conservative.f90
DEP_F90_CONSE=\
	".\Debug\my_messages.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\fluids.f90
DEP_F90_FLUID=\
	".\Debug\Riemann_barotropic_Phi.mod"\
	".\Debug\Riemann_polytropic_Phi.mod"\
	".\Debug\Riemann_van_der_waals_Phi.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\main.f90
NODEP_F90_MAIN_=\
	".\Debug\Riemann_solver.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\my_messages.f90
# End Source File
# Begin Source File

SOURCE=.\polytropic.f90
DEP_F90_POLYT=\
	".\Debug\conserved_quantities.mod"\
	".\Debug\my_messages.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Riemann_barotropic.f90
DEP_F90_RIEMA=\
	".\Debug\barotropic.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Riemann_polytropic.f90
DEP_F90_RIEMAN=\
	".\Debug\polytropic.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Riemann_solver.f90
DEP_F90_RIEMANN=\
	".\Debug\fluids.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Riemann_van_der_waals.f90
DEP_F90_RIEMANN_=\
	".\Debug\van_der_Waals.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\van_der_waals.f90
DEP_F90_VAN_D=\
	".\Debug\conserved_quantities.mod"\
	".\Debug\my_messages.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
