# Miscellaneous useful functions

export printover


"Copy a type, modifying specified fields"
function cpmod{T}(pp::T, di)
    di = !isa(di, Associative) ? Dict(di) : di
    ns = fieldnames(pp)
    args = Array(Any, length(ns))
    for (i,n) in enumerate(ns)
        args[i] = get(di, n, getfield(pp, n))
    end
    T(args...)
end
cpmod{T}(pp::T; kws...) = cpmod(pp, kws)
# TODO: Remove uses of cpmod once default constructors for immutables are
# changed

"Map a function across rows of a 2D array"
maprows(f, m::Matrix) = mapslices(f, m, 2)

"Is this value not NaN?"
notnan(x) = !isnan(x)

"Does this contain any NaN values?"
hasnan(x) = any(isnan(x))

"Wait for ENTER to be pressed before continuing"
function wait_for_enter()
    println("Press enter↵ to continue...")
    readline(STDIN)
end

"Print over a previous line"
function printover(line)
    if isdefined(Main, :IJulia)
        print("\r" * line)
    else
        print("\u1b[1G")   # go to first column
        print(line)
        print("\u1b[K")    # clear the rest of the line
    end
end
