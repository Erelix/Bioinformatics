with Ada.Text_IO; use Ada.Text_IO;
with Ada.Strings.Unbounded; use Ada.Strings.Unbounded;

package FASTA is

    function Load_FASTA (Input : File_Type) return Unbounded_String;

    function Load_FASTA (File_Name : String) return Unbounded_String;

end FASTA;
