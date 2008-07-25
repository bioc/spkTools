## retrieve the spikein matrix from a SpikeInExpressionSet object
setMethod("spikeIn", "SpikeInExpressionSet",
          function(object)
          assayDataElement(object, "spikeIn")
          )

