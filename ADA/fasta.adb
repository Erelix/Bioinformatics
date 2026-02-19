package body FASTA is

    function Load_FASTA (Input : File_Type) return Unbounded_String is
       Header : String := Get_Line (Input);
       Result : Unbounded_String := Null_Unbounded_String;
    begin
       while not End_Of_File (Input) loop
          declare
             Line : String := Get_Line (Input);
          begin
             exit when Line (1) = '>';
             Append (Result, Line);
          end;
       end loop;
       return Result;
    end;

    function Load_FASTA (File_Name : String) return Unbounded_String is
       Input : File_Type;
       Result : Unbounded_String;
    begin
       Open (Input, In_File, File_Name);
       Result := Load_FASTA (Input);
       Close (Input);
       return Result;
    end;

end;
