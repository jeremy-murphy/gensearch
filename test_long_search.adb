with Text_Io; use Text_Io; 
with Ada.Integer_Text_Io; use Ada.Integer_Text_Io; 
procedure Test_Long_Search is
  F: Integer;
  Number_Of_Tests: Integer;
  Pattern_Size: Integer;
  Increment: Integer;
  type Character_Sequence is array(Integer range <>) of Character;
  type Integer_Sequence is array(Integer range <>) of Integer;
  type Skip_Sequence is array(Character range <>) of Integer;
  
  Max_Size: constant Integer := 200_000;
  C: Character;
  S1: Character_Sequence(0 .. Max_Size);
  Base_Line, I, S1_Length, S2_Length: Integer;
  File: Text_Io.File_Type;
  
  procedure Compute_Next(pattern: Character_Sequence; a, m: Integer; 
                           next: out Integer_Sequence) is
    j: Integer := a;
    t: Integer := a - 1;
  begin
    next(a) := a - 1;
    while j < m - 1 loop
      while t >= a and then pattern(j) /= pattern(t) loop
         t := next(t);
      end loop;
      j := j + 1; t := t + 1;
      if pattern(j) = pattern(t) then
        next(j) := next(t);
      else
        next(j) := t;
      end if;
    end loop;
  end Compute_Next;
  
  function KMP(text, pattern: Character_Sequence; 
               b, n, a, m: Integer) return Integer is
    pattern_size, j, k: Integer;
    next: Integer_Sequence(a .. m - 1);
  begin
    Compute_Next(pattern, a, m, next);
    
    pattern_size := m - a; j := a; k := b;
    while j < m and then k < n loop 
      while j >= a and then text(k) /= pattern(j) loop
        j := next(j); 
      end loop; 
      k := k + 1; j := j + 1; 
    end loop; 
    if j = m then 
      return k - pattern_size; 
    else 
      return n; 
    end if;
    
  end KMP;
  
  function L(text, pattern: Character_Sequence; 
                b, n, a, m: Integer) return Integer is
    pattern_size, j, k: Integer;
    next: Integer_Sequence(a .. m - 1);
  begin
    pattern_size := m - a;
    Compute_Next(pattern, a, m, next);
    
    pattern_size := m - a; k := b;
    if pattern_size = 1 then 
      while k /= n and then text(k) /= pattern(a) loop
        k := k + 1;
      end loop;
      return k;
    end if;
    
    while k <= n - pattern_size loop
      while text(k) /= pattern(a) loop
        k := k + 1; 
        if k > n - pattern_size then 
          return n;
        end if;
      end loop;
      
      j := a + 1; k := k + 1;
      while text(k) = pattern(j) loop
        k := k + 1; j := j + 1; 
        if j = m then 
          return k - pattern_size;
        end if;
      end loop;
      
      loop
        j := next(j);
        if j < a then 
           k := k + 1; exit; 
        end if;
        exit when j = a;
        while text(k) = pattern(j) loop
          k := k + 1; j := j + 1; 
          if j = m then 
            return k - pattern_size;
          end if;
          if k = n then 
            return n;
          end if;
        end loop;
      end loop;
      
    end loop;
    return n;
    
  end L;
  
  function SF(text, pattern: Character_Sequence; 
               b, n, a, m: Integer) return Integer is
    pattern_size, j, k: Integer;
  begin
    pattern_size := m - a; k := b;
    if pattern_size = 1 then 
      while k /= n and then text(k) /= pattern(a) loop
        k := k + 1;
      end loop;
      return k;
    end if;
    
    while k <= n - pattern_size loop
      while text(k) /= pattern(a) loop
        k := k + 1; 
        if k > n - pattern_size then 
          return n;
        end if;
      end loop;
      
      j := a + 1; k := k + 1;
      while text(k) = pattern(j) loop
        k := k + 1; j := j + 1; 
        if j = m then 
          return k - pattern_size;
        end if;
      end loop;
      
      k := k - (j - a) + 1;
    end loop;
    return n;
  end SF;
  
  function AL(text, pattern: Character_Sequence; 
              b, n, a, m: Integer) return Integer is
    pattern_size, text_size, j, k, large, adjustment, mismatch_shift: Integer;
    next: Integer_Sequence(a .. m - 1);
    skip: Skip_Sequence(Character'Range);
  begin
    pattern_size := m - a; text_size := n - b; k := b;
    if pattern_size = 1 then 
      while k /= n and then text(k) /= pattern(a) loop
        k := k + 1;
      end loop;
      return k;
    end if;
    
    Compute_Next(pattern, a, m, next);
    
    for i in Character'Range loop
      skip(i) := pattern_size;
    end loop;
    for j in a .. m - 2 loop
      skip(pattern(j)) := m - 1 - j;
    end loop;
    mismatch_shift := skip(pattern(m - 1));
    skip(pattern(m - 1)) := 0;
    
    large := text_size + 1;
    skip(pattern(m - 1)) := large;
    adjustment := large + pattern_size - 1;
    k := k - n;
    loop
      k := k + pattern_size - 1;
      exit when k >= 0;
      loop
        k := k + skip(text(n + k)); 
        exit when k >= 0;
      end loop;
      if k < pattern_size then
        return n;
      end if;
      k := k - adjustment;
      
      if text(n + k) /= pattern(a) then
        k := k + mismatch_shift;
      else 
        j := a + 1;
        loop
          k := k + 1;
          exit when text(n + k) /= pattern(j);
          j := j + 1; 
          if j = m then 
            return n + k - pattern_size + 1;
          end if;
        end loop;
        
        if mismatch_shift > j - a then
          k := k + (mismatch_shift - (j - a));
        else
          loop
            j := next(j);
            if j < a then 
               k := k + 1; 
               exit; 
            end if;
            exit when j = a;
            while text(n + k) = pattern(j) loop
              k := k + 1; j := j + 1; 
              if j = m then 
                return n + k - pattern_size;
              end if;
              if k = 0 then 
                return n;
              end if;
            end loop;
          end loop;
          
        end if;
      end if;
      
    end loop;
    return n;
    
  end AL;
  
  subtype hash_range is Integer range 0..255;
  
  function hash(text: Character_Sequence; k: Integer) return hash_range;
  pragma inline(hash);
  
  function hash(text: Character_Sequence; k: Integer) return hash_range is
  begin
    return hash_range(character'pos(text(k)));
  end hash;
  
  suffix_size: constant Integer := 1;
  
  function HAL(text, pattern: Character_Sequence; 
                 b, n, a, m: Integer) return Integer is
    pattern_size, text_size, j, k, large, adjustment, mismatch_shift: Integer;
    next: Integer_Sequence(a .. m - 1);
    skip: Integer_Sequence(hash_range);
  begin
    pattern_size := m - a; text_size := n - b; k := b;
    if pattern_size = 1 then 
      while k /= n and then text(k) /= pattern(a) loop
        k := k + 1;
      end loop;
      return k;
    end if;
    
    Compute_Next(pattern, a, m, next);
    
    
    for i in hash_range loop
      skip(i) := pattern_size - suffix_size + 1;
    end loop;
    for j in a + suffix_size - 1 .. m - 2 loop
      skip(hash(pattern, j)) := m - 1 - j;
    end loop;
    mismatch_shift := skip(hash(pattern, m - 1));
    skip(hash(pattern, m - 1)) := 0;
    
    large := text_size + 1;
    skip(hash(pattern, m - 1)) := large;
    adjustment := large + pattern_size - 1;
    k := k - n;
    loop
      k := k + pattern_size - 1;
      exit when k >= 0;
      loop
        k := k + skip(hash(text, n + k)); 
        exit when k >= 0;
      end loop;
      if k < pattern_size then
        return n;
      end if;
      k := k - adjustment;
      
      if text(n + k) /= pattern(a) then
        k := k + mismatch_shift;
      else 
        j := a + 1;
        loop
          k := k + 1;
          exit when text(n + k) /= pattern(j);
          j := j + 1; 
          if j = m then 
            return n + k - pattern_size + 1;
          end if;
        end loop;
        
        if mismatch_shift > j - a then
          k := k + (mismatch_shift - (j - a));
        else
          loop
            j := next(j);
            if j < a then 
               k := k + 1; 
               exit; 
            end if;
            exit when j = a;
            while text(n + k) = pattern(j) loop
              k := k + 1; j := j + 1; 
              if j = m then 
                return n + k - pattern_size;
              end if;
              if k = 0 then 
                return n;
              end if;
            end loop;
          end loop;
          
        end if;
      end if;
      
    end loop;
    return n;
    
  end HAL;
  
  
  function AL0(text, pattern: Character_Sequence; 
                 b, n, a, m: Integer) return Integer is
    pattern_size, j, k, d, mismatch_shift: Integer;
    next: Integer_Sequence(a .. m - 1);
    skip: Skip_Sequence(Character'Range);
  begin
    pattern_size := m - a; k := b;
    if pattern_size = 1 then 
      while k /= n and then text(k) /= pattern(a) loop
        k := k + 1;
      end loop;
      return k;
    end if;
    
    Compute_Next(pattern, a, m, next);
    
    for i in Character'Range loop
      skip(i) := pattern_size;
    end loop;
    for j in a .. m - 2 loop
      skip(pattern(j)) := m - 1 - j;
    end loop;
    mismatch_shift := skip(pattern(m - 1));
    skip(pattern(m - 1)) := 0;
    
    while k <= n - pattern_size loop
      loop
        d := skip(text(k + pattern_size - 1)); 
        exit when d = 0;
        k := k + d; 
        if k > n - pattern_size then 
          return n;
        end if;
      end loop;
      
      j := a; 
      while text(k) = pattern(j) loop
        k := k + 1; j := j + 1; 
        if j = m - 1 then 
          return k - pattern_size + 1;
        end if;
      end loop;
      
      if mismatch_shift > j - a then
        k := k + (mismatch_shift - (j - a));
      else
        loop
          j := next(j);
          if j < a then 
             k := k + 1; exit; 
          end if;
          exit when j = a;
          while text(k) = pattern(j) loop
            k := k + 1; j := j + 1; 
            if j = m then 
              return k - pattern_size;
            end if;
            if k = n then 
              return n;
            end if;
          end loop;
        end loop;
        
      end if;
    end loop;
    return n;
    
  end AL0;
  
  function SF1(text, pattern: Character_Sequence; 
               b, n, a, m: Integer) return Integer is
    pattern_size, j, k, k0: Integer;
  begin
    pattern_size := m - a;
    if n < m then
      return n;
    end if;
    j := a; k := b; k0 := k;
    while j /= m loop
      if text(k) /= pattern(j) then
        if k = n - pattern_size then
          return n;
        else
          k0 := k0 + 1; k := k0; j := a;
        end if;
      else
        k := k + 1; j := j + 1;
      end if;
    end loop;
    return k0;
  end SF1;
  
  function SF2(text, pattern: Character_Sequence; 
               b, n, a, m: Integer) return Integer is
    pattern_size, j, k, k0, n0: Integer;
  begin
    pattern_size := m - a;
    if n - b < pattern_size then
      return n;
    end if;
    j := a; k := b; k0 := k; n0 := n - b;
    while j /= m loop
      if text(k) = pattern(j) then
        k := k + 1; j := j + 1;
      else
        if n0 = pattern_size then
          return n;
        else
          k0 := k0 + 1; k := k0; j := a; n0 := n0 - 1;
        end if;
      end if;
    end loop;
    return k0;
  end SF2;
  
  type Algorithm_Enumeration is (Dummy, SF, SF1, SF2, L, AL, HAL);
  
    Algorithm_Names: array(Algorithm_Enumeration) of String(1 .. 17) :=
      ("selection code   ", 
       "SF               ", 
       "HP SF            ", 
       "SGI SF           ", 
       "L                ",
       "AL               ",
       "HAL              ");
    
  function Algorithm(k: Algorithm_Enumeration; 
                     text, pattern: Character_Sequence;
                     b, n, a, m: Integer) return Integer is
  begin
    case k is
      when Dummy => return b;
      when SF => return SF(text, pattern, b, n, a, m);
      when SF1 => return SF1(text, pattern, b, n, a, m);
      when SF2 => return SF2(text, pattern, b, n, a, m);
      when L => return L(text, pattern, b, n, a, m);
      when AL => return AL(text, pattern, b, n, a, m);
      when HAL => return HAL(text, pattern, b, n, a, m);
    end case;
  end Algorithm;
  
  procedure Report(K: Algorithm_Enumeration;
                   S1, S2: Character_Sequence; b, n, a, m: Integer) is
    P: Integer;
  begin
    P := Algorithm(K, S1, S2, b, n, a, m);
    Put("    String "); Put('"');
    for I in a .. m - 1 loop
      Put(S2(I));
    end loop;
    
    if P = n then
      Put(" not found");
      New_Line;
    else
      Put('"'); Put(" found at position ");
      Put(P);
      New_Line;
    end if;
    if Base_Line = 0 then
      Base_Line := P - b;
    else
      if P - b /= Base_Line then
        Put("*****Incorrect result!"); New_Line;
      end if;
    end if;
  end Report;
  
  S2: Character_Sequence(0 .. 100);
begin
  Put("Input Number of tests and pattern size: "); Text_Io.Flush;
  Get(Number_Of_Tests); 
  Get(Pattern_Size);
  New_Line; Put("Number of tests: "); Put(Number_Of_Tests); New_Line;
  Put("Pattern size: "); Put(Pattern_Size); New_Line;
  S2_Length := Pattern_Size;
  
  Text_Io.Open(File, Text_IO.In_File, "long.txt");
  Text_Io.Set_Input(File);
  
  I := 0;
  while not Text_Io.End_Of_File loop
    Text_Io.Get_Immediate(C);
    S1(I) := C;
    I := I + 1;
  end loop;
  S1_Length := I;
  Put(S1_Length); Put(" characters read."); New_Line;
  
  Increment := (S1_Length - S2_Length) / Number_Of_Tests;
  F := 0;
  for K in 1 .. Number_Of_Tests loop
    for I in 0 .. Pattern_Size - 1 loop
      S2(I) := S1(F + I);
    end loop;
    F := F + Increment;
    
    Base_Line := 0;
    for K in Algorithm_Enumeration'Succ(Algorithm_Enumeration'First) .. 
             Algorithm_Enumeration'Last loop
      Put("  Using "); Put(Algorithm_Names(k)); New_Line;
      Report(K, S1, S2, 0, S1_Length, 0, S2_Length);
    end loop;
    New_Line;
    
  end loop;
  
end Test_Long_Search;
