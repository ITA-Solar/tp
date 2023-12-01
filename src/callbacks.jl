
abstract type AbstractTPCallback end

struct DiscreteTPCallback{F1, F2} <: AbstractTPCallback
    condition::F1
    affect!  ::F2
    function DiscreteTPCallback(condition::F1, affect!::F2) where {F1, F2}
        new{F1, F2}(condition, affect!)
    end
end
function DiscreteTPCallback()
    condition(x) = false
    affect! = nothing
    return DiscreteTPCallback(condition, affect!)
end


