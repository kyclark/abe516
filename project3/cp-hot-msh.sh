while read -r ID NAME; do
    echo $ID $NAME
    cp mash/$ID.msh hot/$NAME.msh
done < hot-ids

echo Done.
