function axgxdata = import_axgx(filename)
% function axgxdata = import_axgx(filename)
%%
%% Function imports axgx datafiles (see the bottom of the file for complete
%% doumentation cut and paste from AXOGRAPH_READEWRITE_.h
%% 
%% Input arguments: 
%%      filename: string containing filename to load
%% Output:
%%      axgxdata: struct containing
%%              title: column titles
%%              data:  column data
%%
%%
%%  Bassam Atallah 092706   Only deals with datatypes 9 and 10


fid = fopen(filename, 'r', 'b');
fseek(fid, 0, 'bof')

%% NOTE: alot of space could be saved if data was saved as 'single' rather
%% then the default 'double'

filetype = fread(fid,4,'*char') %%FILE EXTENSION
if ~strcmp(filetype','axgx')
    error('Axograph import only compatible with axgx files')
end
fread(fid,1,'long') %% VERSION

NCol = fread(fid,1,'long') %% NCOL (NSWEEPS)
for i=1:NCol
    np = fread(fid,1,'long') %% Npoints (totalnumber of colums)
    datatype= fread(fid,1,'long') %% Number of bytes in String
    title_length = fread(fid,1,'long') %%
    axgxdata(i).title= fread(fid,title_length,'*char') %%

    switch(datatype)
        case 9 %% series
            value = fread(fid,1,'double') %%
            incre = fread(fid,1,'double') %%
            axgxdata(i).data = [value:incre:np*(incre-1)];
        case 10 %% scaled short
            scal = fread(fid,1,'double') %%
            off = fread(fid,1,'double') %%
            axgxdata(i).data = fread(fid,np,'short') %%
            axgxdata(i).data = axgxdata(i).data*scal+ off;
        otherwise
            error('AXOGRAPH file of Unknown Datatype')
    end

end
%% check if end of file
comment_length = fread(fid,1,'long') %%
comments= fread(fid,comment_length,'*char') %%
note_length = fread(fid,1,'long') %%
note= fread(fid,note_length,'*char') %%

fclose(fid);


% The three AxoGraph file formats are descibed here for completeness, but this  
% information is not needed in order to use the supplied AxoGraph file functions
% to read and write binary data. For information about reading and writing graph
% display information, see the end of the section on the AxoGraph X file format. 
% 
% 
% AxoGraph Data File Format
% =========================
% 
% Header
% ------
% Byte	Type		Contents
% 0		char[4]		AxoGraph file header identifier = 'AxGr' - same as document type ID
% 4		short		AxoGraph graph file format ID = 1 
% 6		short		Number of columns to follow
% 
% 
% Each column
% -----------
% Byte	Type		Contents
% 0		long		Number of points in the column ( columnPoints )
% 4		char[80]	Column title (Pascal 'String' format) - S.I. units should be in brackets e.g. 'Current (pA)'
% 84		float		1st Data point
% 88		float		2nd Data point
% ..		..			....
% ..		..			etc.
% 
% 
% ----------------------------------------------------------------------------------
% 
% AxoGraph Digitized Data File Format
% ===================================
% 
% Header
% ------
% Byte	Type		Contents
% 0		char[4]		AxoGraph file header identifier = 'AxGr' - same as document type ID
% 4		short		AxoGraph file format ID = 2
% 6		short		Number of columns to follow
% 
% 
% Each column
% ----------------------
% Byte	Type		Contents
% 0		long		Number of points in the column ( columnPoints )
% 4		long		Data type
% 8		char[80]	Column title (Pascal 'String' format) - S.I. units should be in brackets e.g. 'Current (pA)'
% 84		float		Scaling Factor
% 88		short		1st Data point
% 90		short		2nd Data point
% ..		...			....
% ..		...			etc.
% 
% 
% ----------------------------------------------------------------------------------
% 
% AxoGraph X Data File Format
% ===================================
% 
% Header
% ------
% Byte	Type		Contents
% 0		char[4]		AxoGraph file header identifier = 'AxGx' - same as filename extension
% 4		long		AxoGraph X file format ID = 3
% 8		long		Number of columns to follow
% 
% 
% Each column
% ----------------------
% Byte	Type		Contents
% 0		long		Number of points in the column ( columnPoints )
% 4		long		Column type 
% 8		long		Length of column title in bytes (Unicode - 2 bytes per character)
% 12		char*		Column title (Unicode 2 byte per char) - S.I. units should be in brackets e.g. 'Current (pA)'
% ??		??			Byte offset depends on length of column title string. 
% ..		...			Numeric type and layout depend on the column type
% ..		...			....
% ..		...			etc.
% 
% 
% Six column types are supported...
% 	4: short 
% 	5: long
% 	6: float
% 	7: double
% 	9: 'series'
% 	10: 'scaled short'
% 
% In the first four column types, data is stored as a simple array of the corresponding type.
% The 'scaled short' column type stores data as a 'double' scaling factor and offset, and a 'short' array.
% The 'series' column type stores data as a 'double' first value and a 'double' increment.
% 
% Prior to AxoGraph X, all graph display information was stored in the 'resource fork' of the file, 
% and the resource fork format was not documented. In contrast, AxoGraph X has a 'flat' format
% with all display information stored immediately following the data columns. 
% It is safe to simply leave out this information. AxoGraph X will use default parameters 
% when the file is read in. For greater control of graph appearance when creating a file 
% it may be necessary to add display format information. When reading in a file,  
% it may be necessary to access the 'Notes' string. The following is a preliminary description 
% of the file format used to store important elements of graph display information. 
% It is not supported in the AxoGraph_ReadWrite example functions. 
% 
% The Comment and Notes strings are stored immediately after the last data column.
% Both are stored in Unicode string format..
% 
% 
% Unicode string format
% ----------------------
% 	long		Length of string in bytes
% 	char*		Notes string (Unicode 2 byte per char)
% 				For Latin1 strings, every second byte is an ASCII character code
% 
% Each trace consists of a pair of columns. The trace header specifies the 
% X and Y column numbers, and other trace-specific information. 
% 'bool' header fields are stored as long int: false = 0, true = 1
% 
% The number of traces is stored immediately after the comment and notes strings.
% 
% 	long		Number of trace headers to follow
% 
% Header for each trace
% ----------------------
% 	long		X column number
% 	long		Y column number
% 	long		Error bar column number or -1 if no error bars
% 
% 	long		Group number that this column belongs to
% 	bool		Trace shown? False if trace is hidden
% 
% 	double		Minimum X data point in this trace
% 	double		Maximum X data point in this trace (if both are zero, they will be recalculated)
% 	double		Minimum positive X data point in this trace (used in log-axis format)
% 
% 	bool		True if X axis data is regularly spaced
% 	bool		True if X axis data is monotonic (each point > previous point)
% 	double		Interval between points for regular X axis data
% 
% 	double		Minimum Y data point in this trace
% 	double		Maximum Y data point in this trace (if both are zero, they will be recalculated)
% 	double		Minimum positive Y data point in this trace (used in log-axis format)
% 
% 	long		Trace color with RGB values serialized into a long int
% 
% 	bool		True if a line plot joining the data points is displayed
% 	double		Thickness of the line plot (can be less than 1.0 for fine lines)
% 	long		Pen style (zero for solid line, non zero for dashed lines)
% 
% 	bool		True if symbols are displayed
% 	long		Symbol type
% 	long		Symbol size (radius in pixels)
% 
% 	bool		True if some symbols are to be skipped 
% 	bool		True if symbols are to be skipped by distance instead of number of points
% 	long		Minimum separation of symbols in pixes is previous parameter is true
% 
% 	bool		True for a histogram plot
% 	long		Type of histogram (zero for standard solid fill)
% 	long		Separation between adjacent histogram bars expressed as a percentage of bar width
% 
% 	bool		True if error bars are displayed 
% 	bool		True if a positive error bar is displayed 
% 	bool		True if a negative error bar is displayed 
% 	long		Error bar width in pixels
% 
% ---------------------------------------------------------------------------------- */


