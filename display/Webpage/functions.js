

    //since javascript doesn't have a built-in range function	
function range(start, stop){
	var myList = [];
	for (let i = start; i <= (stop-1); i++){
		myList.push(i);
	}
	
    return myList;
    }



    //Operator JavaScript


    let letterColors = {"A":"red", "T":"green", "C":"blue","G":"black"};
    var operator_container = document.getElementById("operator_container")	


    function renderOperator(operator_data) {
    
        let subcontainer = document.createElement("div");
        subcontainer.classList.add("operator_subcontainer");

        let info_container = document.createElement("div");
        info_container.classList.add("operator_info_container");
        
        let perc_ident = document.createElement("div");
        perc_ident.classList.add("operator_info");
        perc_ident.innerText = "% identity cutoff: "+operator_data["perc_ident"];
        info_container.appendChild(perc_ident);

        let input_seq = document.createElement("div");
        input_seq.classList.add("operator_info");
        input_seq.innerText = "Input sequence: "+operator_data["input_seq"];
        info_container.appendChild(input_seq);
        
        let num_seqs = document.createElement("div");
        num_seqs.classList.add("operator_info");
        num_seqs.innerText = "Sequences compared: "+operator_data["num_seqs"];
        info_container.appendChild(num_seqs);

        let operator = document.createElement("div");
        operator.classList.add("operator");

        let letters = operator_data["motif"].map(x => x["base"]);
        let scale = operator_data["motif"].map(x => Math.pow(x["score"],3));		

        for (i in letters) {
            let letter = document.createElement('div');
            letter.innerText = letters[i];
            letter.setAttribute("style", "color:"+letterColors[letters[i]]+"; font-size:4rem; transform:scaleY("+scale[i].toString()+"); display:inline-block; padding-bottom: 2.5rem; width:3rem;")
            operator.appendChild(letter);
        }
        subcontainer.appendChild(info_container);	
        subcontainer.appendChild(operator);	
        operator_container.appendChild(subcontainer);
    }




	// Operon Javascript


    let all_operons = document.getElementById("all_operons_container");
	
    let template_container = document.createElement("div");
    template_container.classList.add("operonContainer");

    let homolog_container = document.createElement("div");
    homolog_container.classList.add("operon_scroll_container");



	function renderOperon(meta_data, operon) {
 
    let operon_container = document.createElement("div");

    let operon_graphic = document.createElement("div");
    operon_graphic.classList.add("operon_graphic");
    let operon_description = document.createElement("div");
    operon_description.classList.add("operon_description");
    //operon_description.innerText = "Accession: "+accession+"\n\nPercent Identity: "+identity;
    operon_description.innerText = "Accession: "+meta_data["accession"]+"\n\nPercent Identity: "+meta_data["identity"]+"\n\nOrganism: "+meta_data["organism"];

	let container = document.createElement("div");
    container.classList.add("operonContainer");

	for (i in operon) {
		
		let link = document.createElement('a');
		if (operon[i]['link'].length !== 0){
			link.setAttribute("href", "https://ncbi.nlm.nih.gov/protein/"+operon[i]['link']); }
		link.setAttribute("target", "_blank");
		let gene = document.createElement('div');
		if (operon[i]["direction"] == "+"){
			gene.classList.add("Rgene");
		} else {
			gene.classList.add("Lgene");
		}

		gene.setAttribute("style", "background-color: "+operon[i]["color"]+
			"; width:"+operon[i]["length"]+
			"; margin-right:"+operon[i]["spacer"])
		
		gene.style.setProperty('--color', '20px solid '+operon[i]["color"]);
		
		description = document.createElement('span');
		description.innerText = operon[i]["description"];
		description.classList.add("description");
		gene.appendChild(description)
		link.appendChild(gene);

		operon_graphic.appendChild(link);

	}
    let spacer = document.createElement('span');
	spacer.classList.add("spacer");

	operon_container.appendChild(operon_description);
	operon_container.appendChild(operon_graphic);
	operon_container.appendChild(spacer);
    
    if (meta_data["identity"] == "100.0"){
    template_container.appendChild(operon_container);
    all_operons.appendChild(template_container);
    } else {
    homolog_container.appendChild(operon_container);
    all_operons.appendChild(homolog_container);
	}

    }
