
-------------------------------------------------------
-- Lua functions for annotating VCFs using vcfanno
-- Version number is linked to the version of the toml
-------------------------------------------------------

-- Splits 'str' into pieces matching 'pat', returns them as an array
-- from: http://lua-users.org/wiki/SplitJoin
local function split(str,pat)
   
    local tbl = {}
    str:gsub(pat, function(x) tbl[#tbl+1]=x end)
    return tbl
end

-- Applies a function to given tbl with numeric keys 
-- modifies values in place
local function apply_operation(operation, tbl)

    for k, v in ipairs(tbl) do
        tbl[k] = operation(v)
    end
end

-- Takes annotation from SpliceAI and returns max DS
-- Annotation format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
function spliceai_max_score(spliceai_anno)

    local pattern = string.format("([^%s]+)", '|')
    local split_fields = split(spliceai_anno[1], pattern )
    local score_table = { unpack(split_fields, 3, 6) }
    apply_operation(tonumber, score_table)

    return math.max(unpack(score_table))
end
